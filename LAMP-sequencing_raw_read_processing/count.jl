files=ARGS

using CodecZlib
using Statistics

SEQS = ("AAGATGGTATTTCT", "AGCTGGACTTCCCT", "ATGGTGCTAACAA")
COLS = ("Plate", "Well", "UMI", "Seq1", "Seq2", "Seq3", "Count")
COLS_KMERS = ("Plate", "Well", "Kmer", "Count")
KMER_K = 31
KMER_DIFF = 11

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

function newPath(path::AbstractString, suffix::AbstractString, ext = nothing, dir = nothing)
    (path, ex) = greedysplitext(path)
    ext = isnothing(ext) ? ex : ext
    (pathdir, file) = splitdir(path)
    dir = isnothing(dir) ? pathdir : dir
    dir = dir == "" ? "./" : dir
    if !isdir(dir)
        # throw(ErrorException("OUTDIR does not exist."))
        println(string("'", dir, "' does not exist. Creating."))
        mkpath(dir)
    end
    return string(dir, "/", file, suffix, ext)
end

function format_count(x::Int64, total::Int64)
    return "$x\t($(round(x/total*1000)/10) %)"
end

function retrieve_key(readname::AbstractString)
    (rname, rdes) = split(readname, ' ')
    rn = split(rname, ':')
    return rn[end-2:end]
end

function match_SEQS(readseq::AbstractString)
    ms = [findfirst(s, readseq) for s in SEQS]
    ms = [!isnothing(m) ? m[1] : 0 for m in ms]
    join(ms, ':')
end

function extract_kmers(seq::AbstractString, k::Int64)
    n = length(seq) - k + 1
    kmers = Vector{AbstractString}(undef, n)
    for i in 1:n
        kmers[i] = seq[i:(i+k-1)]
    end
    return unique(kmers)
end

function update_kmer_count(kmers, counts, minmaxdiff)
    for kmer in kmers
        if !haskey(  counts, kmer )
            counts[ kmer ] = 1
        else
            counts[ kmer ] += 1
        end
    end
    ## clean-up counts
    mincounts = minimum(values(counts))
    maxcounts = maximum(values(counts))
    for k in keys(counts)
        if counts[k] <= maxcounts - minmaxdiff
            delete!(counts, k)
        end
    end
    return counts
end

for file in files
    counts = Dict{String, Int64}()
    kmercounts = Dict{ String, Dict{String, Int64} }()
    stats = Dict("TOTAL" => 0,
                 "FAIL" => 0)

    println( file )
	  f = GzipDecompressorStream( open( file ) )
	  while !eof(f)
	      r = [ readline(f) for i in 1:4 ]
        rk = retrieve_key( r[1] )
        rs = match_SEQS( r[2] )
        site = string( join(rk, ':'), ':', rs)
        if !haskey( counts, site )
            counts[ site ] = 1
        else
            counts[ site ] += 1
        end
        if rs == "0:0:0"
            kmers = extract_kmers( r[2], KMER_K )
            kmersite = join( rk[1:2], ':' )
            if !haskey( kmercounts, kmersite )
                kmercounts[ kmersite ] = Dict(k => 1 for k in kmers)
            else
                kmercounts[ kmersite ] = update_kmer_count(kmers, kmercounts[ kmersite ], KMER_DIFF )
            end
            stats["FAIL"] += 1
        end
        stats["TOTAL"] += 1
        if mod( stats["TOTAL"], 1000000 ) == 0
            println(stats["TOTAL"])
	      end
    end
    f2 = open(newPath(file, "_counted", ".tsv", "Counted"), "w")
    println(f2, join(COLS, '\t'))
    ## sort counts and write
    for k in sort(collect(keys(counts)))
        println(f2, join(vcat(split(k, ':'), counts[k]), '\t'))
    end

    f3 = open(newPath(file, "_unmatchedkmers", ".tsv", "Counted"), "w")
    println(f3, join(COLS_KMERS, '\t'))
    for k in sort(collect(keys( kmercounts )))
        for (kmer, count) in sort(collect( kmercounts[k] ), by=x->x[2], rev=true)
            println(f3, join(vcat(split(k, ':'), kmer, count), '\t'))
        end
    end

    ## report final statistics
    println(stderr, string(stats["TOTAL"], "\treads were read."))
    println(stderr, string(format_count(stats["FAIL"], stats["TOTAL"]), " reads were not counted."))

    close(f)
    close(f2)
    close(f3)
end
