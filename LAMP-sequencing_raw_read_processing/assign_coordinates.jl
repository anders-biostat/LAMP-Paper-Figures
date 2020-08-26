using ArgParse
using BioSequences
using CodecZlib
using Statistics

function parse_commandline()
    s = ArgParseSettings("Simple Demux script handling offline UMIs\n" *
                         "Konrad Herbst, k.herbst@zmbh.uni-heidelberg.de," *
                         "Sophie Herbst",
                         version = "Version 0.1",
                         add_version = true)

    @add_arg_table! s begin
        "--n_max_errors", "-e"
            help = "maximal errors allowed for demultiplexing"
            arg_type = Int
            default = 0
        "--distance", "-d"
            help  = "distance metric used for demultiplexing" *
                "(allowed values: hamming OR seqlev)"
            range_tester = (x->x=="hamming"||x=="seqlev")
            default = "hamming"
            arg_type = String
        "--outdir", "-o"
            help = "directory to write files into"
            default = ""
            arg_type = String
        "indices_P7"
            help = "Indices of P7 site: Plate identity (csv format: name, sample, seq_bc)"
            required = true
            arg_type = String
        "indices_P5"
            help = "Indices of P5 site: Well identity (csv format: name, sample, seq_bc)"
            required = true
            arg_type = String
        "readfiles"
            help = "read1 fastq/fastq.gz file(s)"
            required = true
            nargs = 'R'
            arg_type = String
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

# ## TODO for interactive testing
# parsed_args = Dict(
#     "n_max_errors" => 2,
#     "distance" => :hamming,
#     "indices_P7" => "indices-P7-plates.csv",
#     "indices_P5" => "indices-P5-wells.csv",
#     "readfiles" => ["head100k-Undetermined_S0_R1_001.fastq.gz"],
#     "outdir" => ""
# )

if parsed_args["distance"] == "hamming"
    parsed_args["distance"] = :hamming
elseif parsed_args["distance"] == "seqlev"
    parsed_args["distance"] = :levenshtein
else
    error("Wrong distance argument")## should not happen
end

## helper function
function greedysplitext(path::AbstractString)
    ext = ""
    next = "."
    while next != ""
        (path, next) = splitext(path)
        ext = string(next, ext)
    end
    return (path, ext)
end

function newPath(path::AbstractString, suffix::AbstractString, dir = nothing)
    (path, ext) = greedysplitext(path)
    (pathdir, file) = splitdir(path)
    dir = dir == nothing ? pathdir : dir
    if !isdir(dir)
        # throw(ErrorException("OUTDIR does not exist."))
        println(string("'", dir, "' does not exist. Creating."))
        mkpath(dir)
    end
    return string(dir, "/", file, suffix, ext)
end

function load_bcs( filename::String )
    samples = Vector{String}()
    seqs = Vector{LongSequence{DNAAlphabet{4}}}()
    for l = eachline( filename )
        f = split( l, "," )
        if l == "name,sample,seq_bc"
            continue
        end
        (name, sample, seq) = split(l, ',')
        push!(samples, sample)
        push!(seqs, LongDNASeq(seq))
    end
    return (samples, seqs)
end

function format_count(x::Int64, total::Int64)
    return "$x\t($(round(x/total*1000)/10) %)"
end

## load BCs and create demux
(samples_P7, seqs_P7) = load_bcs(parsed_args["indices_P7"])
dplxr_P7 = Demultiplexer(seqs_P7, n_max_errors = parsed_args["n_max_errors"], distance = parsed_args["distance"])
(samples_P5, seqs_P5) = load_bcs(parsed_args["indices_P5"])
dplxr_P5 = Demultiplexer(seqs_P5, n_max_errors = parsed_args["n_max_errors"], distance = parsed_args["distance"])

for readfile in parsed_args["readfiles"]
    stats = Dict("TOTAL" => 0,
                 "ASSIGNED" => 0,
                 "FAIL_P7_MISS" => 0,
                 "FAIL_P5_MISS" => 0)
    println(stderr, readfile)
    ## initialize file handles
    if parsed_args["outdir"] == ""
        outdir = string(splitdir(readfile)[1], ".")
    else
        outdir = parsed_args["outdir"]
    end
    if splitext(readfile)[2] == ".gz"
        reads = GzipDecompressorStream(open(readfile, "r"))
        reads_pass = GzipCompressorStream(open(newPath(readfile, "_pass", outdir), "w"))
        reads_fail = GzipCompressorStream(open(newPath(readfile, "_fail", outdir), "w"))
    else
        reads = open(readfile, "r")
        reads_pass = open(newPath(readfile, "_pass", outdir), "w")
        reads_fail = open(newPath(readfile, "_fail", outdir), "w")
    end

    while !eof(reads)
        r = [ readline(reads) for i in 1:4 ]
        (rname, rdes) = split(r[1], ' ')
        (bc_i7, bc_i5) = split(split(rdes, ':')[4], '+')
        umi = bc_i5[12:20]
        bc_i5 = bc_i5[1:11]
        r_i7 = demultiplex(dplxr_P7, LongDNASeq(bc_i7))
        r_i5 = demultiplex(dplxr_P5, LongDNASeq(bc_i5))
        if (r_i7[2] == -1) && (r_i5[2] == -1)
            stats["FAIL_P7_MISS"] += 1
            stats["FAIL_P5_MISS"] += 1
            write(reads_fail, string(join(r, '\n'), '\n'))
        elseif r_i7[2] == -1
            stats["FAIL_P7_MISS"] += 1
            write(reads_fail, string(join(r, '\n'), '\n'))
        elseif r_i5[2] == -1
            stats["FAIL_P5_MISS"] += 1
            write(reads_fail, string(join(r, '\n'), '\n'))
        else
            rname = string(rname, ':', samples_P7[r_i7[1]], ':', samples_P5[r_i5[1]], ':', umi)
            r[1] = string(rname, ' ', rdes)
            write(reads_pass, string(join(r, '\n'), '\n'))
            stats["ASSIGNED"] += 1
        end
        stats["TOTAL"] += 1

        ## report progress once in a while
        if (stats["TOTAL"] > 0) & (stats["TOTAL"] % 100000 == 0)
            println(stderr, string(stats["TOTAL"], " reads processed."))
        end
    end

    ## report final statistics
    println(stderr, string(stats["TOTAL"], "\treads were read."))
    println(stderr, string(format_count(stats["ASSIGNED"], stats["TOTAL"]), " reads were assigned."))
    println(stderr, string(format_count(stats["FAIL_P7_MISS"], stats["TOTAL"]), " reads which could not be assigned based on i7."))
    println(stderr, string(format_count(stats["FAIL_P5_MISS"], stats["TOTAL"]), " reads which could not be assigned based on i5."))

    ## close file handles
    close(reads)
    close(reads_pass)
    close(reads_fail)
end
