#!/usr/bin/env julia
#
#Use this method to build an index of the positions that can become nullomers due to a single bp mutation. 


using ArgParse, BioSequences, BioSequences.FASTA, Bio.Seq, DelimitedFiles

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
      "--nullomer_length", "-l"
        help = "nullomer length"
        arg_type = Int
        default = 11
      "--chromosome", "-c"
        help = "chromosome number (1:23)"
        arg_type = Int
        default = 0
      "--nullomer_start", "-s"
        help = "start of nullomer"
        arg_type = Int
        default = 1
      "--nullomer_end", "-e"
        help = "end of nullomer"
        arg_type = Int
        default = 0
      "--genome_location", "-g"
        help = "path to genome .fa file"
        default = "/lustre/scratch117/cellgen/team218/lh20/cellranger_genomes/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
      "--nullomer_path", "-p"
        help = "path to directory containing nullomer list"
        default = "/lustre/scratch117/cellgen/team218/MH/nullomers/"
      "--nullomer_file", "-f"
        help = "file containing nullomer list"
        default = ""
    end
    return parse_args(s)
end

parsed_args = parse_commandline();
nullomerLen = parsed_args["nullomer_length"];
chr = parsed_args["chromosome"];
nullomerStart = parsed_args["nullomer_start"];
nullomerEnd = parsed_args["nullomer_end"];
genomeLocation = parsed_args["genome_location"];
nullomerPath = parsed_args["nullomer_path"];
nullomerFile = parsed_args["nullomer_file"];

#read the list of nullomers
if nullomerFile==""
    nullomersStrs = convert(Array{String, 1}, readdlm(nullomerPath * "nullomers" * string(nullomerLen) * ".txt")[:,1]);
else
    nullomersStrs = convert(Array{String, 1}, readdlm(nullomerPath * nullomerFile, comments=true)[:,1]);
end
nullomers = DNAKmer[];
for n in nullomersStrs
    push!(nullomers, DNAKmer(n));
end
if (nullomerStart>1) || (nullomerEnd>0) #only consider a subset of the nullomers
    nullomers = nullomers[nullomerStart:min(length(nullomers), nullomerEnd)];
end
reader = FASTA.Reader(open(genomeLocation, "r"));
savefilename = nullomerPath * "nullomers" * string(nullomerLen) * "tmp/nullomers" * string(nullomerLen) * "mutpos.tsv";
if nullomerFile!=""
    savefilename = nullomerPath * replace(nullomerFile, ".txt" => "_mutpos.tsv");
end
if chr!=0
    savefilename = replace(savefilename, ".tsv" => "_" * string(chr) * ".tsv");
    savefilename = replace(savefilename, "tmp/" => "tmp/" * string(chr) * "/");
    if nullomerStart!=0
        savefilename = replace(savefilename, ".tsv" => "_" * string(nullomerStart) * "_" * string(nullomerEnd) * ".tsv");
    end
end
savefilehandle = open(savefilename, "w");
if (chr<=1) && (nullomerStart<=1)
    println(savefilehandle, "#chr\tpos\tmutation\tkmerpos\tnullomer");
end
for rec in reader
    chrStr = seqname(rec);
    c = tryparse(Int64, chrStr);
    if chrStr=="X"
        c = 23;
    elseif chrStr=="Y"
        c = 24;
    end
    if (c!=nothing) && ((chr==0) || (c==chr))
        savefilehandle = open(savefilename, "a"); 
        chrSeq = sequence(rec);
        #find matches
        seqLen = length(chrSeq);
        ninds = findall(chrSeq.!=DNA_N);
        chrSeqN = DNASequence(string(chrSeq)[ninds]);
        seqLenN = length(chrSeqN);
        for j in 1:length(nullomers)
            nullomer = nullomers[j];
            nullomerQuery = ApproximateSearchQuery(DNASequence(nullomer)); #preprocessing to speed up
            start = ninds[1];
            s = 1;
            while (start>0) && (start<seqLen-nullomerLen) && (s>0) && (s<seqLenN)
                s = approxsearchindex(chrSeqN, nullomerQuery, 1, s);
                if (s<=0) || (s>=seqLenN)
                    break;
                end
                start = ninds[s];
                if composition(chrSeq[start:start+nullomerLen-1])[DNA_N]==0
                    #Find out the exact location of the mismatch and what mutation it was caused by
                    nullomerDNA = DNASequence(convert(String, nullomer));
                    for i in 1:nullomerLen
                        if chrSeq[start-1+i]!=nullomer[i]
                            mutLoc = start + i - 1;
                            #look at the remaining bases to decide if this was an insert, a deletion or a substitution
                            origBase = string(chrSeq[start-2+i:start+i]);
                            newBase = string(nullomerDNA[i])
                            if chrSeq[start+i:start+nullomerLen]==nullomerDNA[i:nullomerLen] #deletion
                                newBase = "-";
                            elseif chrSeq[start+i-1:start+nullomerLen-2]==nullomerDNA[i+1:nullomerLen] #insertion
                                origBase = string(chrSeq[start-2+i]) * "-" * string(chrSeq[start-1+i]);
                                mutLoc -= 1; 
                            end
                            println(savefilehandle, chrStr * "\t" * string(mutLoc) * "\t" * origBase * ":" * newBase * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                            if (i==1) && (string(origBase[2])!=string("-")) && (string(origBase[1])==string(newBase)) #we can also delete the first aa since it is the same as the previous one
                                println(savefilehandle, chrStr * "\t" * string(mutLoc) * "\t" * string(chrSeq[start-2+i:start+i]) * ":-\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                            end
                            if i==nullomerLen #if the insertion or deletion happened at the last base, then it should also work as a substitution. Note that we do not need to do this correction for the first aa since the loop automatically finds those cases
                                if string(origBase[2])==string("-")
                                    println(savefilehandle, chrStr * "\t" * string(mutLoc + 1) * "\t" * string(chrSeq[start-2+i:start+i]) * ":" * newBase * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j)); #remember to compensate for the change in position for the insertion above
                                end
                                if (string(origBase[3])==newBase) #we can also delete the last base since it is the same as the next one
                                    println(savefilehandle, chrStr * "\t" * string(mutLoc + 1) * "\t" * string(chrSeq[start-2+i:start+i]) * ":-\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                elseif string(newBase)==string("-") #we can do a substitution or insertion to achieve the same effect as the deletion
                                    println(savefilehandle, chrStr * "\t" * string(mutLoc) * "\t" * string(chrSeq[start-2+i:start+i]) * ":" * string(chrSeq[start+i]) * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                    println(savefilehandle, chrStr * "\t" * string(mutLoc) * "\t" * string(chrSeq[start-2+i]) * "-" * string(chrSeq[start-1+i]) * ":" * string(chrSeq[start+i]) * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                end
                            end
                            break;
                        end
                    end
                end
                s += 1;
            end
        end
        close(savefilehandle);
    end
end
close(reader);

