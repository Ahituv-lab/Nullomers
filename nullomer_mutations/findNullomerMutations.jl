#Use this method to scan the mutations reported in a patient or Gnomad to find out which ones result in nullomers
#

using DelimitedFiles, ArgParse, StatsBase, Plots, Bio.Intervals

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
      "--nullomer_length", "-l"
        help = "nullomer length"
        arg_type = Int
        default = 11
      "--chromosome", "-c"
        help = "chromosome"
        arg_type = Int
        default = 0
      "--gnomad", "-g"
        help = "consider the gnomad variants instead of the mutations"
        arg_type = Bool
        default = false
      "--nullomer_file", "-f"
        help = "file with nullomer locations, relative to nullomer_path"
        arg_type = String
        default = ""
      "--population", "-p"
        help = "which population to consider"
        arg_type = String
        default = ""
      "--nullomer_path", "-P"
        help = "path to directory containing nullomer substitution list"
	default = "/lustre/scratch117/cellgen/team218/MH/nullomers/"
    end
    return parse_args(s)
end

parsed_args = parse_commandline();
nullomerLen = parsed_args["nullomer_length"];
chromosome = parsed_args["chromosome"];
gnomad = parsed_args["gnomad"];
population = parsed_args["population"];
nullomerPath = parsed_args["nullomer_path"];
nullomerfile = parsed_args["nullomer_file"];

chrcol = 1;
poscol = 2;
mutcol = 5;

#load the file containing the positions of the mutations
suffix = "sub";
savefilename = nullomerPath * "nullomers" * string(nullomerLen) * suffix * "s_gnomad.tsv";
if population!=""
    savefilename = nullomerPath * "nullomers" * string(nullomerLen) * suffix * "s_gnomad_" * population * ".tsv";
end
if chromosome>0
    savefilename = replace(savefilename, ".tsv" => "_" * string(chromosome) * ".tsv");
end

mutfilename = "/lustre/scratch117/cellgen/team218/MH/nullomers/vcfs/Gnomad/snps_only_gnomad_afraw.vcf";
if population!=""
    mutfilename = "/lustre/scratch117/cellgen/team218/MH/nullomers/vcfs/gnomad_" * population * ".vcf";
end
if (chromosome>0) && (population=="")
    mutfilename = replace(mutfilename, ".vcf" => "_" * string(chromosome) * ".vcf");
end

@time muts = readdlm(mutfilename, comments=true);
muts[:,chrcol] = [replace(i,"chr"=>"") for i in muts[:,chrcol]];
xinds = findall((muts[:,chrcol].=="X"));
yinds = findall((muts[:,chrcol].=="Y"));
muts[xinds, chrcol] = 23*ones(Int64, length(xinds), 1);
muts[yinds, chrcol] = 24*ones(Int64, length(yinds), 1);
somaticinds = setdiff(1:size(muts,1),union(xinds,yinds));
muts[somaticinds, chrcol] = tryparse.(Int64, muts[somaticinds,chrcol]);
if chromosome>0
    muts = muts[findall(muts[:,chrcol].==chromosome),:];
end

nullomerfilename = nullomerPath * "nullomers" * string(nullomerLen) * suffix * "spos.tsv";
if (nullomerLen>12) && (population=="") && gnomad
    nullomerfilename = nullomerPath * "nullomers13tmp/nullomers" * string(nullomerLen) * suffix * "spos_" * string(chromosome) * ".tsv"
end
if nullomerfile!=""
    nullomerfilename = nullomerPath * nullomerfile;
    savefilename = nullomerPath * replace(nullomerfile, ".tsv" => "_gnomad.tsv");
    if chromosome>0
        savefilename = replace(savefilename, ".tsv" => "_" * string(chromosome) * ".tsv");
    end
end
@time nullomermuts = readdlm(nullomerfilename, comments=true);
xinds = findall((nullomermuts[:,1].=="X"));
nullomermuts[xinds, 1] = 23*ones(Int64, length(xinds), 1);
yinds = findall((nullomermuts[:,1].=="Y"));
nullomermuts[yinds, 1] = 24*ones(Int64, length(yinds), 1);

pops = ["fin", "nfe", "asj", "oth", "eas", "sas", "afr", "ami", "amr"];
println("Saving results to " * savefilename)
savefilehandle = open(savefilename, "w");
if chromosome<=1
    if population==""
        println(savefilehandle, "#chr\tpos\tmutation\tvariant\taf\tnullomer");
    else
        print(savefilehandle, "#chr\tpos\tmutation\tvariant");
        for p in pops
            print(savefilehandle, "\taf_" * p);
        end
        println(savefilehandle, "\tnullomer");
    end
end

for i in 1:24
    #find all of the mutations that result in nullomers
    mutChrInds = findall(muts[:,chrcol].==i);
    nullomerChrInds = findall(nullomermuts[:,1].==i);
    tmp = indexin(nullomermuts[nullomerChrInds,2], muts[mutChrInds,poscol]);
    nullomerMutInds = findall(tmp.!=nothing);
    tmp2 = indexin(muts[mutChrInds,poscol], nullomermuts[nullomerChrInds,2]);
    mutInds = findall(tmp2.!=nothing);
    if length(mutInds)>0
        if population==""
            mutintervalcollection = IntervalCollection{Tuple{String, String, Float64}}([Interval("chr" * string(muts[mutChrInds[mutInds[j]],chrcol]), muts[mutChrInds[mutInds[j]],poscol], muts[mutChrInds[mutInds[j]],poscol], '.', (string(muts[mutChrInds[mutInds[j]],3]), string(muts[mutChrInds[mutInds[j]],5]), muts[mutChrInds[mutInds[j]],6])) for j in 1:length(mutInds)], true);
        else
            mutintervalcollection = IntervalCollection{Tuple{String, String, Float64, String}}([Interval("chr" * string(muts[mutChrInds[mutInds[j]],chrcol]), muts[mutChrInds[mutInds[j]],poscol], muts[mutChrInds[mutInds[j]],poscol], '.', (string(muts[mutChrInds[mutInds[j]],3]), string(muts[mutChrInds[mutInds[j]],5]), muts[mutChrInds[mutInds[j]],6], string(muts[mutChrInds[mutInds[j]],8]))) for j in 1:length(mutInds)], true);
        end
        nullomerintervalcollection = IntervalCollection{Tuple{String, Int64}}([Interval("chr" * string(nullomermuts[nullomerChrInds[nullomerMutInds[j]],1]), nullomermuts[nullomerChrInds[nullomerMutInds[j]],2], nullomermuts[nullomerChrInds[nullomerMutInds[j]],2], '.', (string(nullomermuts[nullomerChrInds[nullomerMutInds[j]],3]), nullomermuts[nullomerChrInds[nullomerMutInds[j]],5])) for j in 1:length(nullomerMutInds)], true)
        #found locations that can result in nullomers, now also make sure that we have the right type of mutation
        for it in eachoverlap(mutintervalcollection, nullomerintervalcollection)
            if population!=""
                tokens = split(metadata(it[1])[4], ";");
                #extract the information about the allele frequencies
                afs = Float64[];
                for p in pops
                    afstrind = findfirst(startswith.(tokens, "AF_" * p * "="));
                    af = nothing;
                    if afstrind!=nothing
                        af = tryparse(Float64, convert(Array{String, 1}, split(split(tokens[afstrind], "=")[2], ","))[1]);
                    end
                    if af==nothing
                        af = 0.0;
                    end
                    push!(afs, af);
                end
                if string(metadata(it[2])[1][end])==string(metadata(it[1])[2])
                    println(savefilehandle, string(i) * "\t" * string(leftposition(it[1])) * "\t" * metadata(it[2])[1] * "\t" * metadata(it[1])[1] * "\t" * join(string.(round.(afs; digits=6)), "\t") * "\t" * string(metadata(it[2])[2]));
                end
            else
                if string(metadata(it[2])[1][end])==string(metadata(it[1])[2])
                    println(savefilehandle, string(i) * "\t" * string(leftposition(it[1])) * "\t" * metadata(it[2])[1] * "\t" * metadata(it[1])[1] * "\t" * string(round(metadata(it[1])[3]; digits=6)) * "\t" * string(metadata(it[2])[2]));
                end
            end
        end
    end
end
close(savefilehandle);


