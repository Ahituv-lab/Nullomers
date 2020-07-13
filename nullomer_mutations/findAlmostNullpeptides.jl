#!/usr/bin/env julia
#
#Use this method to build an index of the positions that can become nullpeptides due to a single bp mutation. 
#The problem with this map is that it considers all possible AA changes. If we only want to consider the ones that can occur through single bp substitution, then we need to do more work. 


using ArgParse, BioSequences, BioSequences.FASTA, Bio.Seq, DelimitedFiles

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
      "--nullomer_length", "-l"
        help = "nullomer length"
        arg_type = Int
        default = 4
      "--nullomer_start", "-s"
        help = "start of nullomer"
        arg_type = Int
        default = 1
      "--nullomer_end", "-e"
        help = "end of nullomer"
        arg_type = Int
        default = 0
      "--proteome_location", "-p"
        help = "path to proteome .fa file"
        default = "/path/to/file/in/fasta/format"
      "--nullomer_path", "-n"
        help = "path to directory containing nullomer list"
        default = "/path/to/directory/"
    end
    return parse_args(s)
end

parsed_args = parse_commandline();
nullomerLen = parsed_args["nullomer_length"];
nullomerStart = parsed_args["nullomer_start"];
nullomerEnd = parsed_args["nullomer_end"];
proteomeLocation = parsed_args["proteome_location"];
nullomerPath = parsed_args["nullomer_path"];

#read the list of nullomers
nullpeptideStrs = convert(Array{String, 1}, readdlm(nullomerPath * "nullpeptides" * string(nullomerLen) * ".txt")[:,1]);
nullpeptides = AminoAcidSequence[];
for n in nullpeptideStrs
    push!(nullpeptides, AminoAcidSequence(n));
end
if (nullomerStart>1) || (nullomerEnd>0)
    nullpeptides = nullpeptides[nullomerStart:min(length(nullpeptides), nullomerEnd)];
end
reader = FASTA.Reader(open(proteomeLocation, "r"));
suffix = "subs";
savefilename = nullomerPath * "nullpeptides" * string(nullomerLen) * suffix * "pos.tsv";
if nullomerEnd!=0
    savefilename = nullomerPath * "nullpeptides" * string(nullomerLen) * "tmp/nullpeptides" * string(nullomerLen) * suffix * "pos_" * string(nullomerStart) * "_" * string(nullomerEnd) * ".tsv";
end
savefilehandle = open(savefilename, "w");
if (nullomerStart<=1)
    println(savefilehandle, "#transcriptid\tpos\tmutation\taapos\tnullpeptide");
end
for rec in reader
    tokens = split(identifier(rec), "|");
    ind = findfirst(startswith.(tokens, "ENST"));
    if (ind!=nothing) && ((transcriptid==nothing) || findfirst(startswith.(tokens, transcriptid))!=nothing)
        genename = tokens[ind];
        aas = FASTA.sequence(rec);
        aaLen = length(aas);
        for j in 1:length(nullpeptides)
            nullpeptide = nullpeptides[j];
            start = 1;
            while (start>0) && (start<aaLen-nullomerLen)
                start = approxsearchindex(aas, nullpeptide, 1, start);
                if (start<=0) || (start>=aaLen-nullomerLen)
                    break;
                end
                #Find out the exact location of the mismatch and what mutation it was caused by
                for i in 1:nullomerLen
                    if (start>1) && (aas[start-1+i]!=nullpeptide[i])
                        mutLoc = start + i - 1;
                        #look at the remaining bases to decide if this was an insert, a deletion or a substitution
                        origAA = string(aas[start-2+i:start+i]);
                        newAA = string(nullpeptide[i])
                        if aas[start+i:start+nullomerLen]==nullpeptide[i:nullomerLen] #deletion
                            newAA = "-";
                        elseif aas[start+i-1:start+nullomerLen-2]==nullpeptide[i+1:nullomerLen] #insertion
                            origAA = string(aas[start-2+i]) * "-" * string(aas[start-1+i]);
                            mutLoc -= 1;
                        end
                        if indels || (!indels && (newAA!=string("-")) && (string(origAA[2])!=string("-")))
                            println(savefilehandle, genename * "\t" * string(mutLoc) * "\t" * origAA * ":" * newAA * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                        end
                        if indels && (i==1) && (string(origAA[2])!=string("-")) && (string(origAA[1])==string(newAA)) #we can also delete the first aa since it is the same as the previous one
                            println(savefilehandle, genename * "\t" * string(mutLoc) * "\t" * string(aas[start-2+i:start+i]) * ":-\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                        end
                        if i==nullomerLen #if the insertion or deletion happened at the last base, then it should also work as a substitution. Note that we do not need to do this correction for the first aa since the loop automatically finds those cases
                            if string(origAA[2])==string("-")
                                println(savefilehandle, genename * "\t" * string(mutLoc + 1) * "\t" * string(aas[start-2+i:start+i]) * ":" * newAA * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j)); #remember to compensate for the change in position for the insertion above
                                if indels && (string(origAA[3])==newAA) #we can also delete the last base since it is the same as the next one
                                    println(savefilehandle, genename * "\t" * string(mutLoc + 1) * "\t" * string(aas[start-2+i:start+i]) * ":-\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                end
                            elseif string(newAA)==string("-") #we can do a substitution or insertion to achieve the same effect as the deletion
                                println(savefilehandle, genename * "\t" * string(mutLoc) * "\t" * string(aas[start-2+i:start+i]) * ":" * string(aas[start+i]) * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                if indels
                                    println(savefilehandle, genename * "\t" * string(mutLoc) * "\t" * string(aas[start-2+i]) * "-" * string(aas[start-1+i]) * ":" * string(aas[start+i]) * "\t" * string(i) * "\t" * string(nullomerStart - 1 + j));
                                end
                            end
                        end
                        break;
                    end
                end
                start += 1;
            end
        end
    end
end
close(savefilehandle);
close(reader);

