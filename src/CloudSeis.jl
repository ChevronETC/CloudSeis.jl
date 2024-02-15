module CloudSeis

using AbstractStorage, Blosc, Distributed, JSON, Pkg, Random, TeaSeis, UUIDs

struct TracePropertyDef{T}
    label::String
    description::String
    format::Type{T}
    elementcount::Int
end
TracePropertyDef(p::Dict) = TracePropertyDef(p["label"], p["description"], stringtype2type(p["format"], p["elementcount"]), p["elementcount"])

struct TraceProperty{T}
    def::TracePropertyDef{T}
    byteoffset::Int
end
TraceProperty(p::Dict) = TraceProperty(TracePropertyDef(p), p["byteoffset"])

Base.Dict(tracepropdef::TracePropertyDef) = Dict(
    "label"=>tracepropdef.label,
    "description"=>tracepropdef.description,
    "format"=>propertyformatstring(tracepropdef),
    "elementcount"=>tracepropdef.elementcount)

Base.Dict(traceprop::TraceProperty) = merge(Dict(traceprop.def), Dict("byteoffset"=>traceprop.byteoffset))

struct DataProperty{T}
    label::String
    value::T
end

Base.Dict(dataproperty::DataProperty) = Dict(dataproperty.label=>dataproperty.value)
DataProperty(dataproperty::Dict) = DataProperty(first(keys(dataproperty)), first(values(dataproperty)))

include("stockprops.jl")

struct Geometry
    u1::Int
    un::Int
    v1::Int
    vn::Int
    w1::Int
    wn::Int
    ox::Float64
    oy::Float64
    oz::Float64
    ux::Float64
    uy::Float64
    uz::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    wx::Float64
    wy::Float64
    wz::Float64
end

function Geometry(;
        ox=0.0,oy=0.0,oz=0.0,
        ux=1.0,uy=0.0,uz=0.0,
        vx=0.0,vy=1.0,vz=0.0,
        wx=0.0,wy=0.0,wz=1.0,
        u1=1,un=2,
        v1=1,vn=2,
        w1=1,wn=2)
    Geometry(u1,un,v1,vn,w1,wn,ox,oy,oz,ux,uy,uz,vx,vy,vz,wx,wy,wz)
end

Base.Dict(g::Geometry) = Dict(
    "ox"=>g.ox, "oy"=>g.oy, "oz"=>g.oz,
    "ux"=>g.ux, "uy"=>g.uy, "uz"=>g.uz,
    "vx"=>g.vx, "vy"=>g.vy, "vz"=>g.vz,
    "wx"=>g.wx, "wy"=>g.wy, "wz"=>g.wz,
    "u1"=>g.u1, "un"=>g.un,
    "v1"=>g.v1, "vn"=>g.vn,
    "w1"=>g.w1, "wn"=>g.wn)

struct Extent{C<:Container}
    name::String
    container::C
    frameindices::UnitRange{Int}
end

function Dict(extent::Extent)
    Dict("name"=>extent.name, "container"=>minimaldict(extent.container), "firstframe"=>extent.frameindices[1], "lastframe"=>extent.frameindices[end])
end

function Extent(extent::Dict, containers::Vector{C}) where {C<:Container}
    if haskey(extent, "container")
        if !haskey(extent["container"], "prefix") # for backwards compatability
            extent["container"]["prefix"] = ""
        end
        Extent(extent["name"], Container(C, extent["container"], session(containers[1])), extent["firstframe"]:extent["lastframe"])
    else # for backwards compatability -- reading data generated with CloudSeis <= 0.2
        Extent(extent["name"], containers[1], extent["firstframe"]:extent["lastframe"])
    end
end

Base.length(extent::Extent) = length(extent.frameindices)

const CACHE_NONE = 0
const CACHE_FOLDMAP = 1
const CACHE_ALL = 2
const CACHE_ALL_LEFT_JUSTIFY = 3

# compression algorithms <--
abstract type AbstractCompressor end

function pkgversion(uuid)
    pkgs = Pkg.dependencies()
    pkgs[uuid].version
end

struct BloscCompressor <: AbstractCompressor
    algorithm::String
    level::Int
    shuffle::Bool
end
BloscCompressor(d::Dict) = BloscCompressor(d["library_options"]["algorithm"], d["library_options"]["level"], d["library_options"]["shuffle"])
function Dict(c::BloscCompressor)
    Dict(
        "method" => "blosc",
        "library" => "Blosc.jl",
        "library_version" => string(pkgversion(UUID("a74b3585-a348-5f62-a45c-50e91977d574"))), # this UUID is specific to Blosc.jl
        "library_options" => Dict(
            "algorithm" => c.algorithm,
            "level" => c.level,
            "shuffle" => c.shuffle))
end
Base.copy(c::BloscCompressor) = BloscCompressor(c.algorithm, c.level, c.shuffle)
hdrlength_multipleof(_::BloscCompressor) = 1

struct ZfpCompressor{T,P,R} <: AbstractCompressor
    tol::T
    precision::P
    rate::R
end
function ZfpCompressor(;kwargs...)
    if length(kwargs) == 0
        return ZfpCompressor(nothing, nothing, nothing)
    end

    length(kwargs) == 1 || error("zfp options are at most one of ':tol', ':precision', and ':rate'")
    key,value = keys(kwargs)[1],values(kwargs)[1]

    if key == :tol
        return ZfpCompressor(value, nothing, nothing)
    elseif key == :precision
        return ZfpCompressor(nothing, value, nothing)
    elseif key == :rate
        return ZfpCompressor(nothing, nothing, value)
    else
        error("zfp options are at most one of ':tol', ':precision', and ':rate'")
    end
end

kwargs(compressor::ZfpCompressor{<:Real,Nothing,Nothing}) = (tol=compressor.tol,)
kwargs(compressor::ZfpCompressor{Nothing,<:Real,Nothing}) = (precision=compressor.precision,)
kwargs(compressor::ZfpCompressor{Nothing,Nothing,<:Real}) = (rate=compressor.rate,)
kwargs(_::ZfpCompressor{Nothing,Nothing,Nothing}) = ()

function ZfpCompressor(d::Dict)
    if haskey(d["library_options"], "tol")
        return ZfpCompressor(Float64(d["library_options"]["tol"]), nothing, nothing)
    elseif haskey(d["library_options"], "precision")
        return ZfpCompressor(nothing, Int32(d["library_options"]["precision"]), nothing)
    elseif haskey(d["library_options"], "rate")
        return ZfpCompressor(nothing, nothing, Int(d["library_options"]["rate"]))
    else
        return ZfpCompressor(nothing, nothing, nothing)
    end
end

function Dict(c::ZfpCompressor)
    local library_options
    if c.tol !== nothing
        library_options = Dict("tol" => c.tol)
    elseif c.precision !== nothing
        library_options = Dict("precision" => c.precision)
    elseif c.rate !== nothing
        library_options = Dict("rate" => c.rate)
    else
        library_options = Dict()
    end

    Dict(
        "method" => "zfp",
        "library" => "ZfpCompression.jl",
        "library_version" => string(pkgversion(UUID("43441a71-1662-41c6-b8ea-40ed1525242b"))), # this UUID is specific to ZfpCompression.jl
        "library_options" => library_options
    )
end
Base.copy(c::ZfpCompressor) = ZfpCompressor(c.tol, c.precision, c.rate)
hdrlength_multipleof(_::ZfpCompressor) = 4

struct CloudSeisCvxCompressor <: AbstractCompressor
    b1::Int
    b2::Int
    b3::Int
    scale::Float32
end
function CloudSeisCvxCompressor(;b1=16,b2=16,b3=16,scale=0.001f0)
    if b1 > 1
        b1 = clamp(nextpow(2, b1), 8, 256)
    end
    if b2 > 1
        b2 = clamp(nextpow(2, b2), 8, 256)
    end
    if b3 > 1
        b3 = clamp(nextpow(2, b3), 8, 256)
    end
    CloudSeisCvxCompressor(b1,b2,b3,scale)
end
kwargs(c::CloudSeisCvxCompressor) = (b1=c.b1, b2=c.b2, b3=c.b3, scale=scale)
CloudSeisCvxCompressor(d::Dict) = CloudSeisCvxCompressor(d["library_options"]["b1"], d["library_options"]["b2"], d["library_options"]["b3"], d["library_options"]["scale"])
function Dict(c::CloudSeisCvxCompressor)
    Dict(
        "method" => "cvx",
        "library" => "CvxCompress.jl",
        "library_version" => string(pkgversion(UUID("a74b3585-a348-5f62-a45c-50e91977d574"))), # this UUID is specific to CvxCompress.jl
        "library_options" => Dict(
            "b1" => c.b1,
            "b2" => c.b2,
            "b3" => c.b3,
            "scale" => c.scale
        )
    )
end

Base.copy(c::CloudSeisCvxCompressor) = CloudSeisCvxCompressor(c.b1, c.b2, c.b3, c.scale)
hdrlength_multipleof(_::CloudSeisCvxCompressor) = 4

struct LeftJustifyCompressor <: AbstractCompressor end
LeftJustifyCompressor(d::Dict) = LeftJustifyCompressor()
Dict(c::LeftJustifyCompressor) = Dict("method" => "leftjustify")
Base.copy(c::LeftJustifyCompressor) = LeftJustifyCompressor()
hdrlength_multipleof(_::LeftJustifyCompressor) = 1

struct NotACompressor <: AbstractCompressor end
NotACompressor(d::Dict) = NotACompressor()
Dict(c::NotACompressor) = Dict("method" => "none")
Base.copy(c::NotACompressor) = NotACompressor()
hdrlength_multipleof(_::NotACompressor) = 1

function Compressor(d::Dict)
    method = get(d, "method", "") # enables backwards compatability
    if method == "blosc"
        return BloscCompressor(d)
    elseif method == "zfp"
        return ZfpCompressor(d)
    elseif method == "cvx"
        return CloudSeisCvxCompressor(d)
    elseif method == "leftjustify"
        return LeftJustifyCompressor(d)
    else
        return NotACompressor(d)
    end
end
# -->

mutable struct Cache{Z<:AbstractCompressor}
    extentindex::Int
    data::Vector{UInt8}
    compressor::Z
    type::Int
end
Cache(compressor::AbstractCompressor) = Cache(0, UInt8[], compressor, CACHE_NONE)
Base.similar(cache::Cache) = Cache(0, UInt8[], copy(cache.compressor), CACHE_NONE)

struct CSeis{T,N,Z<:AbstractCompressor,C<:Container,U<:NamedTuple,V<:NamedTuple,W<:NamedTuple}
    containers::Vector{C}
    mode::String
    datatype::String
    traceformat::Type{T}
    byteorder::String
    axis_units::NTuple{N,String}
    axis_domains::NTuple{N,String}
    axis_lengths::NTuple{N,Int}
    axis_pstarts::NTuple{N,Float64}
    axis_pincs::NTuple{N,Float64}
    axis_lstarts::NTuple{N,Int}
    axis_lincs::NTuple{N,Int}
    axis_propdefs::U
    traceproperties::V
    dataproperties::W
    geometry::Union{Geometry,Nothing}
    extents::Vector{Extent{C}}
    cache::Cache{Z}
    hdrlength::Int
end

function Base.copy(io::CSeis, mode, extents)
    CSeis(
        io.containers,
        mode,
        io.datatype,
        io.traceformat,
        io.byteorder,
        io.axis_units,
        io.axis_domains,
        io.axis_lengths,
        io.axis_pstarts,
        io.axis_pincs,
        io.axis_lstarts,
        io.axis_lincs,
        io.axis_propdefs,
        io.traceproperties,
        io.dataproperties,
        io.geometry,
        extents,
        similar(io.cache),
        io.hdrlength)
end

"""
    io = csopen(containers[, mode="r"]; axis_lengths::Vector{Int}[, optional keyword arguments])

returns a handle to a CloudSeis data-set, and where `container::Container`
corresponds to a dataset where all extents are in a single container, and 
`container::Vector{<:Container}` corresponds to a dataset where extents
are sharded across mulitple containers.  The container(s) can be either POSIX
folders (container::Folder), or cloud storage container (e.g. container::AzContainer).

`mode` is one of `"r"` - read, `"w"` - write, `"r+"` - open existing data-set for reading and writing.

If `mode="w"`, then `axis_lenghts` is a required key-word argument.  `axis_lengths` specifies the size of the
container.  The `axis_lengths` vector is of at-least length 3.

# Optional keyword arguments
* `similarto::String` An existing CloudSeis dataset.  If set, then all other named arguments can be used to modify the data context that belongs to the existing CloudSeis dataset.
* `datatype::String` Examples are `CMP`, `SHOT`, etc.  If not set, then `UNKNOWN` is used.
* `traceformat::DataType=Float32` Not supported.  We have only tested against `traceformat=Float32`.
* `byteorder::String`
* `extents::Vector{UnitRange}` List of integer ranges (frame range) for each extent.  Each extent maps to an object in block storage (set only one of `extents`, `frames_per_extent` and `mbytes_per_extent`).
* `frames_per_extent::Int` Nominal number of frames per extent (set only one of `extents`, `frames_per_extent` and `mbytes_per_extent`).
* `mbytes_per_extent::Int` Nominal size of each extent in units of mega-bytes (set only one of `extents`, `frames_per_extent` and `mbytes_per_extent`).
* `geometry::Geometry` An optional three point geometry can be embedded in the CloudSeis file.
* `tracepropertydefs::Vector{TracePropertyDef}` An array of trace property definitions.   These are in addition to the axes property definitions.
* `dataproperties::Vector{DataProperty}` An array of data properties.  One property per dataset rather than one property per trace (as is true for `tracepropertydefs`).
* `axis_propdefs::Vector{TracePropertyDef}` Trace properties corresponding to the CloudSeis axes.  If not set, then `SAMPLE`, `TRACE`, `FRAME`, `VOLUME` and `HYPRCUBE` are used.
* `axis_units::Vector{String}` Units corresponding to the CloudSeis axes. e.g. `SECONDS`, `METERS`, etc.  If not set, then `UNKNOWN` is used.
* `axis_domains::Vector{String}` Domains corresponding to the CloudSeis axes. e.g. `SPACE`, `TIME`, etc.  If not set, then `UNKNOWN` is used.
* `axis_pstarts::Vector{Float64}` Physical origins for each axis.  If not set, then `0.0` is used for the physical origin of each axis.
* `axis_pincs::Vector{Float64}` Physical deltas for each axis.  If not set, then `1.0` is used for the physical delta for each axis.
* `axis_lstarts::Vector{Int}` Logical starts for each axis.  If not set, then `1` is used for the logical start for each axis.
* `axis_lincs::Vector{Int}` Logical increments for each axis.  If not set, then `1` is used for the logical increment for each axis.
* `compressor="leftjustify"` Compress the cache before writing to disk.  This is particularly useful for data with variable fold.  chooose from: ("none", "blosc", "leftjustify", "zfp", "cvx")
* `compressor_options=()` Pass options as key-word arguments to the compression algorithm.  Currently, only `"zfp"` and `"cvx"` have options[1]

# Example
## Azure storage
```julia
using AzStorage, CloudSeis
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container, "w"; axis_lengths=[10,11,12], axis_pincs=[0.004,10.0,20.0])
close(io)
```
## POSIX storage
```
using AzStorage, FolderStorage
container = Folder("mydataset-cs")
io = csopen(container, "w"; axis_lengths=[10,11,12], axis_pincs=[0.004,10.0,20.0])
```

# Notes
* When using the `similarto` option, one can change the number of dimensions of the data-set via `axis_lengths`.  If one shrinks the number of dimensions,
then various data-set properties (e.g. `axis_units`) will be truncated.  The truncation can be customized by using appropriate key-word arguments.

* The zfp compression options are `tol`, `precision` and `rate`.  So, for example:
```
using AzStorage, CloudSeis
container = AzContainer("mydataset-cs"; storageaccount="mystorageacccount")
io = csopen(container, "w"; axis_lengths=[10,11,12], compressor="zfp", compressor_options=(tol=1e-4,))
```
Please refer to ZFPCompressor.jl for more information.  If `compressor_options` is not supplied, then
the defaults are `(precision=16,)`.  Also, note that `compressor_options=()` results
in ZFP lossless compression.  ZFP lossless compression will be used for the headers and foldmap regardless
of the choice of `compressor_options`.

* The cvx compression options are `b1`, `b2`, `b3` and `scale`.  So, for example:
```
```
Please refer to CvxCompress.jl for more information.  If `compressor_options` is no supplied, then
the defaults are `(b1=16,b2=16,b3=16,scale=1e-2)`.
"""
function csopen(containers::Vector{<:Container}, mode;
        similarto = "",
        datatype = "",
        force = false,
        traceformat = nothing,
        byteorder = "",
        extents = UnitRange[],
        frames_per_extent = 0,
        mbytes_per_extent = -1,
        geometry = nothing,
        tracepropertydefs = [],
        dataproperties = [],
        axis_propdefs = [],
        axis_units = [],
        axis_domains = [],
        axis_lengths = [],
        axis_pstarts = [],
        axis_pincs = [],
        axis_lstarts = [],
        axis_lincs = [],
        compressor = nothing,
        compressor_options = nothing)
    kwargs = (
        similarto = similarto,
        datatype = datatype,
        force = force,
        traceformat = traceformat,
        byteorder = byteorder,
        extents = extents,
        frames_per_extent = frames_per_extent,
        mbytes_per_extent = mbytes_per_extent,
        geometry = geometry,
        tracepropertydefs = tracepropertydefs,
        dataproperties = dataproperties,
        axis_propdefs = axis_propdefs,
        axis_units = axis_units,
        axis_domains = axis_domains,
        axis_lengths = axis_lengths,
        axis_pstarts = axis_pstarts,
        axis_pincs = axis_pincs,
        axis_lstarts = axis_lstarts,
        axis_lincs = axis_lincs,
        compressor = compressor,
        compressor_options = compressor_options)
    if mode == "r" || mode == "r+"
        return csopen_read(containers, mode)
    elseif mode == "w" && similarto == ""
        return csopen_write(containers, mode; process_kwargs(;kwargs...)...)
    elseif mode == "w" && similarto != ""
        return csopen_write(containers, mode; process_kwargs_similarto(;kwargs...)...)
    end
    error("mode \"$mode\" is not recognized.")
end
csopen(container::Container, mode; kwargs...) = csopen([container], mode; kwargs...)
csopen(containers::Union{Container,Vector{<:Container}}; kwargs...) = csopen(containers, "r"; kwargs...)

"""
    cscreate(containers[; optional keyword arguments...])

Create a new CloudSeis dataset without creating a corresponding handle and without
opening the data-set.  Please see help for `csopen` for the `optional keyword arguments`.

`containers` is either of type `Container` or `Vector{<:Container}`.  In the former case,
all extents are stored in a single container, and in the later case, extents are sharded
accross multiple containers.

In "r" and "r+" mode, and when the existing data-set is sharded over multiple containers,
it is only required to supply the primary container.  The "description.json" object is
in the "primary" container.
"""
cscreate(containers::Union{Container,Vector{<:Container}}; kwargs...) = close(csopen(containers, "w"; kwargs...))

function csopen_read(containers::Vector{<:Container}, mode)
    description = JSON.parse(read(containers[1], "description.json", String))

    traceproperties = get_trace_properties(description)
    csopen_from_description(containers, mode, description, traceproperties)
end

function make_extents(containers::Vector{C}, nextents, foldername, nominal_frames_per_extent, remaining_frames; last_container_index=0, last_extent_index=0, lastframe=0) where {C<:Container}
    extents = Vector{Extent{C}}(undef, nextents)
    k = last_container_index == length(containers) ? 1 : last_container_index + 1
    l = ceil(Int, log10(last_extent_index + nextents))
    for iextent = 1:nextents
        firstframe = lastframe + 1
        lastframe = firstframe + nominal_frames_per_extent
        if iextent > remaining_frames
            lastframe -= 1
        end
        extents[iextent] = Extent("$foldername/extent-$(lpad(last_extent_index + iextent,l,'0'))", containers[k], firstframe:lastframe)
        k = k == length(containers) ? 1 : k + 1
    end
    extents
end

function csopen_write(containers::Vector{<:Container}, mode; kwargs...)
    ndim = length(kwargs[:axis_lengths])
    ndim == 0 && error("must specify axis_lengths")
    axis_propdefs = kwargs[:axis_propdefs]
    if length(axis_propdefs) == 0
        axis_propdefs = [stockprop[:SAMPLE], stockprop[:TRACE], stockprop[:FRAME], stockprop[:VOLUME], stockprop[:HYPRCUBE]][1:min(ndim,5)]
        for i = 6:ndim
            push!(axis_propdefs, TracePropertyDef("DIM$i", "DIM$i", Int32, 1))
        end
    end
    axis_units = length(kwargs[:axis_units]) == 0 ? ["unknown" for i=1:ndim] : kwargs[:axis_units]
    axis_domains = length(kwargs[:axis_domains]) == 0 ? ["unknown" for i=1:ndim] : kwargs[:axis_domains]
    axis_pstarts = length(kwargs[:axis_pstarts]) == 0 ? [0.0 for i=1:ndim] : kwargs[:axis_pstarts]
    axis_pincs = length(kwargs[:axis_pincs]) == 0 ? [1.0 for i=1:ndim] : kwargs[:axis_pincs]
    axis_lstarts = length(kwargs[:axis_lstarts]) == 0 ? [1 for i=1:ndim] : kwargs[:axis_lstarts]
    axis_lincs = length(kwargs[:axis_lincs]) == 0 ? [1 for i=1:ndim] : kwargs[:axis_lincs]

    traceproperties = get_trace_properties(kwargs[:tracepropertydefs], axis_propdefs)

    traceformat = kwargs[:traceformat]
    axis_lengths = kwargs[:axis_lengths]
    nframes = prod(axis_lengths[3:end])

    local extents
    if length(kwargs[:extents]) > 0
        l = ceil(Int, log10(length(kwargs[:extents])))
        container_index = 0
        nextk() = container_index = container_index == length(containers) ? 1 : container_index+1
        extents = [Extent("extents/extent-$(lpad(iextent,l,'0'))", containers[nextk()], kwargs[:extents][iextent]) for iextent=1:length(kwargs[:extents])]
    else
        local nextents
        if kwargs[:frames_per_extent] > 0
            nextents,r = divrem(prod(axis_lengths[3:end]), kwargs[:frames_per_extent])
            if r > 0
                nextents += 1
            end
        else # mbytes_per_extent
            bytes_total = prod(axis_lengths[2:end])*(axis_lengths[1]*sizeof(traceformat) + headerlength(traceproperties, hdrlength_multipleof(kwargs[:compressor]))) + sizeof(Int)*nframes
            nextents = clamp(div(bytes_total, kwargs[:mbytes_per_extent]*1000*1000), 1, nframes)
        end

        nominal_frames_per_extent, remaining_frames = divrem(nframes, nextents)

        extents = make_extents(containers, nextents, "extents", nominal_frames_per_extent, remaining_frames)
    end

    description = Dict(
        "fileproperties" => Dict(
            "datatype" => kwargs[:datatype],
            "traceformat" => propertyformatstring(traceformat),
            "byteorder" => kwargs[:byteorder],
            "axis_propdefs" => [axis_propdefs[i].label for i=1:ndim],
            "axis_units" => axis_units,
            "axis_domains" => axis_domains,
            "axis_lengths" => axis_lengths,
            "axis_pstarts" => axis_pstarts,
            "axis_pincs" => axis_pincs,
            "axis_lstarts" => axis_lstarts,
            "axis_lincs" => axis_lincs),
        "traceproperties"=>[Dict(traceproperty) for traceproperty in traceproperties],
        "dataproperties"=>isempty(kwargs[:dataproperties]) ? Dict() : mapreduce(Dict, merge, kwargs[:dataproperties]),
        "extents"=>Dict.(extents),
        "compressor"=>Dict(kwargs[:compressor]))

    if kwargs[:geometry] != nothing
        merge!(description, Dict("geometry"=>Dict(kwargs[:geometry])))
    end

    mkpath.(containers)
    write(containers[1], "description.json", json(description, 1))
    io = csopen_from_description(containers, mode, description, traceproperties)
    empty!(io) # agressive emptying (over-writing) of existing data-set with same name.
    io
end

compressors() = Dict("none"=>NotACompressor, ""=>NotACompressor, "blosc"=>()->BloscCompressor("blosclz", 5, true), "zfp"=>ZfpCompressor, "cvx"=>CloudSeisCvxCompressor, "leftjustify"=>LeftJustifyCompressor)
compressor_error() = error("(compressor=$(kwargs[:compressor])) ∉ (\"\",\"none\",\"blosc\",\"zfp\",\"cvx\",\"leftjustify\")")

function process_kwargs(;kwargs...)
    local compressor_options
    if kwargs[:compressor] == "zfp" && kwargs[:compressor_options] === nothing
        compressor_options = (precision=16,)
    elseif kwargs[:compressor_options] === nothing
        compressor_options = ()
    else
        compressor_options = kwargs[:compressor_options]
    end

    compressor = compressors()[kwargs[:compressor] === nothing ? "leftjustify" : kwargs[:compressor]](;compressor_options...)

    (
        similarto = kwargs[:similarto],
        datatype = kwargs[:datatype] == "" ? "custom" : kwargs[:datatype],
        force = kwargs[:force],
        traceformat = kwargs[:traceformat] == nothing ? Float32 : kwargs[:traceformat],
        byteorder = kwargs[:byteorder] == "" ? "LITTLE_ENDIAN" : kwargs[:byteorder],
        extents = kwargs[:extents],
        frames_per_extent = kwargs[:frames_per_extent],
        mbytes_per_extent = kwargs[:mbytes_per_extent] < 0 ? 1024 : kwargs[:mbytes_per_extent],
        geometry = kwargs[:geometry],
        tracepropertydefs = kwargs[:tracepropertydefs],
        dataproperties = kwargs[:dataproperties],
        axis_propdefs = kwargs[:axis_propdefs],
        axis_units = kwargs[:axis_units],
        axis_domains = kwargs[:axis_domains],
        axis_lengths = kwargs[:axis_lengths],
        axis_pstarts = kwargs[:axis_pstarts],
        axis_pincs = kwargs[:axis_pincs],
        axis_lstarts = kwargs[:axis_lstarts],
        axis_lincs = kwargs[:axis_lincs],
        compressor = compressor
    )
end

function process_kwargs_similarto(;kwargs...)
    io = csopen(kwargs[:similarto])

    local mbytes_per_extent, _ndims
    if !(isempty(kwargs[:axis_lengths]))
        @warn "If axis_lengths is set, then the extents will be re-computed.  Use `mbytes_per_extent`, `frames_per_extent`, or `extents` keyword arguments to control the behavior."
        mbytes_per_extent = kwargs[:mbytes_per_extent] < 0 ? 1024 : kwargs[:mbytes_per_extent]
        _ndims = length(kwargs[:axis_lengths])
    else
        mbytes_per_extent = kwargs[:mbytes_per_extent]
        _ndims = ndims(io)
    end

    if _ndims > ndims(io)
        for key in (:axis_propdefs, :axis_units, :axis_domains, :axis_pstarts, :axis_pincs, :axis_lstarts, :axis_lincs)
            if isempty(get(kwargs, key, []))
                error("If increasing the number of dimensions, then all of (axis_propdefs, axis_units, axis_domains, axis_pstarts, axis_pincs, axis_lstarts, axis_lincs) must be specified.")
            end
        end
    elseif _ndims < ndims(io)
        truncating_props = []
        for key in (:axis_propdefs, :axis_units, :axis_domains, :axis_pstarts, :axis_pincs, :axis_lstarts, :axis_lincs)
            if isempty(get(kwargs, key, []))
                push!(truncating_props, key)
            end
        end
        if !isempty(truncating_props)
            @warn "truncating ($(join(truncating_props, ","))) since the number of dimensions is decreased compared to the 'similarto' data-set."
        end
    end

    local extents
    if kwargs[:frames_per_extent] == 0 && mbytes_per_extent < 0 && length(kwargs[:extents]) == 0
        extents = [io.extents[i].frameindices for i=1:length(io.extents)]
    else
        extents = kwargs[:extents]
    end

    names = fieldnames(typeof(io.dataproperties))
    dataproperties = [DataProperty(String(names[i]), io.dataproperties[i]) for i=1:length(io.dataproperties)]

    local compressor
    if kwargs[:compressor] === nothing
        compressor = io.cache.compressor
    else
        local compressor_options
        if kwargs[:compressor] == "zfp" && kwargs[:compressor_options] === nothing
            compressor_options = (precision=16,)
        elseif kwargs[:compressor_options] === nothing
            compressor_options = ()
        else
            compressor_options = kwargs[:compressor_options]
        end
        kwargs[:compressor] ∈ keys(compressors()) || compressor_error()
        compressor = compressors()[kwargs[:compressor]](;compressor_options...)
    end

    (
        similarto = kwargs[:similarto],
        datatype = kwargs[:datatype] == "" ? io.datatype : kwargs[:datatype],
        force = kwargs[:force],
        traceformat = kwargs[:traceformat] == nothing ? io.traceformat : kwargs[:traceformat],
        byteorder = kwargs[:byteorder] == "" ? io.byteorder : kwargs[:byteorder],
        extents = extents,
        frames_per_extent = kwargs[:frames_per_extent],
        mbytes_per_extent = mbytes_per_extent,
        geometry = kwargs[:geometry] == nothing ? io.geometry : kwargs[:geometry],
        tracepropertydefs = isempty(kwargs[:tracepropertydefs]) ? [io.traceproperties[i].def for i=1:length(io.traceproperties)] : kwargs[:tracepropertydefs],
        dataproperties = isempty(kwargs[:dataproperties]) ? dataproperties : kwargs[:dataproperties],
        axis_propdefs = isempty(kwargs[:axis_propdefs]) ? [io.axis_propdefs[i] for i=1:_ndims] : kwargs[:axis_propdefs],
        axis_units = isempty(kwargs[:axis_units]) ? [io.axis_units[i] for i=1:_ndims] : kwargs[:axis_units],
        axis_domains = isempty(kwargs[:axis_domains]) ? [io.axis_domains[i] for i=1:_ndims] : kwargs[:axis_domains],
        axis_lengths = isempty(kwargs[:axis_lengths]) ? [io.axis_lengths[i] for i=1:_ndims] : kwargs[:axis_lengths],
        axis_pstarts = isempty(kwargs[:axis_pstarts]) ? [io.axis_pstarts[i] for i=1:_ndims] : kwargs[:axis_pstarts],
        axis_pincs = isempty(kwargs[:axis_pincs]) ? [io.axis_pincs[i] for i=1:_ndims] : kwargs[:axis_pincs],
        axis_lstarts = isempty(kwargs[:axis_lstarts]) ? [io.axis_lstarts[i] for i=1:_ndims] : kwargs[:axis_lstarts],
        axis_lincs = isempty(kwargs[:axis_lincs]) ? [io.axis_lincs[i] for i=1:_ndims] : kwargs[:axis_lincs],
        compressor = compressor
    )
end

function csopen_from_description(containers, mode, description, traceproperties)
    ndim = length(description["fileproperties"]["axis_lengths"])
    compressor = Compressor(get(description, "compressor", Dict()))
    CSeis(
        containers,
        mode,
        description["fileproperties"]["datatype"],
        stringtype2type(description["fileproperties"]["traceformat"]),
        description["fileproperties"]["byteorder"],
        ntuple(i->description["fileproperties"]["axis_units"][i], ndim),
        ntuple(i->description["fileproperties"]["axis_domains"][i], ndim),
        ntuple(i->description["fileproperties"]["axis_lengths"][i], ndim),
        ntuple(i->description["fileproperties"]["axis_pstarts"][i], ndim),
        ntuple(i->description["fileproperties"]["axis_pincs"][i], ndim),
        haskey(description["fileproperties"], "axis_lstarts") ? ntuple(i->description["fileproperties"]["axis_lstarts"][i], ndim) : ntuple(_->1, ndim), # note backwards compat for data-sets without logical starts/increments
        haskey(description["fileproperties"], "axis_lincs") ? ntuple(i->description["fileproperties"]["axis_lincs"][i], ndim) : ntuple(_->1, ndim), # note backwards compat for data-sets without logical starts/increments
        get_axis_propdefs(description, traceproperties),
        traceproperties,
        get_data_properties(description),
        get_geometry(description),
        [Extent(extent, containers) for extent in description["extents"]],
        Cache(compressor),
        headerlength(traceproperties, hdrlength_multipleof(compressor)))
end

"""
    close(io::CSeis)

Close a handle to a CloudSeis data-set.  This may also flush contents of
the CloudSeis buffer to storage.
"""
function Base.close(io::CSeis)
    if io.mode == "w" || io.mode == "r+"
        flush(io)
    end
    nothing
end

"""
    rm(io::CSeis)

Delete a CloudSeis data-set.
"""
Base.rm(io::CSeis) = rm.(io.containers)

function robust_rm(containers::Vector{<:Container})
    for container in containers
        try
            rm(csopen(container))
        catch
            rm(container)
        end
    end
end

function Base.show(io::IO, cs::CSeis)
    write(io, "CloudSeis file:\n");
    write(io, "\tsize: $(size(cs))\n");
    write(io, "\ttype: $(cs.datatype) ; format: $(propertyformatstring(cs.traceformat))\n");
    write(io, "\taxis domains: $(cs.axis_domains)\n");
    write(io, "\taxis units: $(cs.axis_units)\n");
    write(io, "\taxis properties: $(ntuple(i->cs.axis_propdefs[i].label, length(cs.axis_propdefs)))");
end

function headerlength(traceproperties::NamedTuple, multipleof=1)
    hdrlength = 0
    for traceproperty in traceproperties
        hdrlength += sizeof(eltype(traceproperty.def.format))*traceproperty.def.elementcount
    end
    r = rem(hdrlength, multipleof)
    hdrlength + r
end
"""
    headerlength(io::CSeis)

Returns the length (number of bytes) of a trace header for a CloudSeis data-set.
"""
headerlength(io::CSeis) = io.hdrlength

function get_axis_propdefs(description::Dict, traceproperties)
    c = description["fileproperties"]["axis_propdefs"]
    axis_propdefs = [traceproperties[Symbol(c[i])].def for i=1:length(c)]
    names = ntuple(i->Symbol(axis_propdefs[i].label), length(axis_propdefs))
    NamedTuple{names}(axis_propdefs)
end

function get_trace_properties(description::Dict)
    c = description["traceproperties"]
    traceproperties = [TraceProperty(c[i]) for i=1:length(c)]
    names = ntuple(i->Symbol(traceproperties[i].def.label), length(traceproperties))
    NamedTuple{names}(traceproperties)
end

function get_trace_properties(tracepropdefs, axis_propdefs)
    propdefs = TracePropertyDef[tracepropdefs...]
    for axis_propdef in axis_propdefs
        hasit = false
        for propdef in propdefs
            if axis_propdef.label == propdef.label
                hasit = true
                break
            end
        end
        if !hasit
            push!(propdefs, axis_propdef)
        end
    end

    for std_propdef in [stockprop[:TRC_TYPE]]
        hasit = false
        for propdef in propdefs
            if std_propdef.label == propdef.label
                hasit = true
                break
            end
        end
        if !hasit
            push!(propdefs, std_propdef)
        end
    end

    traceproperties = TraceProperty[]
    byteoffset = 0
    for propdef in propdefs
        push!(traceproperties, TraceProperty(propdef, byteoffset))
        byteoffset += sizeof(eltype(propdef.format))*propdef.elementcount
    end

    names = ntuple(i->Symbol(traceproperties[i].def.label), length(traceproperties))
    NamedTuple{names}(traceproperties)
end

function get_data_properties(description::Dict)
    if "dataproperties" ∉ keys(description)
        return NamedTuple{}()
    end
    c = description["dataproperties"]

    local names, _values
    if isa(c, AbstractArray) # for backwards compatability (datasets generated with version <=1.1.3)
        names = ntuple(i->Symbol(c[i]["label"]), length(c))
        _values = [c[i]["value"] for i=1:length(c)]
    else
        labels = collect(keys(c))
        names = ntuple(i->Symbol(labels[i]), length(labels))
        _values = collect(values(c))
    end
    NamedTuple{names}(_values)
end

function get_geometry(description::Dict)
    if "geometry" ∉ keys(description)
        return nothing
    end

    c = description["geometry"]
    Geometry(
        c["u1"],c["un"],c["v1"],c["vn"],c["w1"],c["wn"],
        c["ox"],c["oy"],c["oz"],c["ux"],c["uy"],c["uz"],
        c["vx"],c["vy"],c["vz"],c["wx"],c["wy"],c["wz"])
end

function stringtype2type(format::AbstractString, elementcount=1)
    format == "INTEGER" && elementcount == 1 && return Int32
    format == "LONG" && elementcount == 1 && return Int64
    format == "FLOAT" && elementcount == 1 && return Float32
    format == "DOUBLE" && elementcount == 1 && return Float64

    format == "INTEGER" && elementcount > 1 && return Vector{Int32}
    format == "LONG" && elementcount > 1 && return Vector{Int64}
    format == "FLOAT" && elementcount > 1 && return Vector{Float32}
    format == "DOUBLE" && elementcount > 1 && return Vector{Float64}
    format == "BYTESTRING" && return Vector{UInt8}
    error("unrecognized format: $(format)")
end

propertyformatstring(T::Union{Type{Int32},Type{Vector{Int32}}}) = "INTEGER"
propertyformatstring(T::Union{Type{Int64},Type{Vector{Int64}}}) = "LONG"
propertyformatstring(T::Union{Type{Float32},Type{Vector{Float32}}}) = "FLOAT"
propertyformatstring(T::Union{Type{Float64},Type{Vector{Float64}}}) = "DOUBLE"
propertyformatstring(T::Type{Vector{UInt8}}) = "BYTESTRING"
propertyformatstring(def::TracePropertyDef) = propertyformatstring(def.format)
propertyformatstring(prop::TraceProperty) = propertyformatstring(prop.def)

function cachesize(io::CSeis, extentindex)
    nframes = length(io.extents[extentindex].frameindices)
    nframes*(sizeof(Int) + (io.hdrlength + sizeof(io.traceformat)*size(io,1))*size(io,2))
end

function cachesize_foldmap(io::CSeis, extentindex)
    nframes = length(io.extents[extentindex].frameindices)
    nframes*sizeof(Int)
end

function notacompressor_cache_from_file!(io::CSeis, extentindex)
    @debug "reading extent $extentindex from block-storage..."
    t = @elapsed begin
        io.cache.data = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end
    mb = length(io.cache.data)/1_000_000
    mbps = mb/t
    @debug "...data read ($mbps MB/s -- $mb MB)"
    nothing
end
function cache_from_file!(io::CSeis{T,N,NotACompressor}, extentindex, regularize) where {T,N}
    regularize || error("'regularize=false' is not supported for 'NotACompressor'")
    notacompressor_cache_from_file!(io, extentindex)
end

function cache_from_file!(io::CSeis{T,N,BloscCompressor}, extentindex, regularize) where {T,N}
    regularize || error("regularize=false is not supported for 'BloscCompressor'")
    @debug "reading and decompressing extent $extentindex..."
    t_read = @elapsed begin
        cdata = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end

    t_decompress = @elapsed begin
        io_cdata = IOBuffer(cdata; read=true, write=false)
        nbuffers = read(io_cdata, Int)
        buffersize = read!(io_cdata, Vector{Int}(undef, nbuffers))

        io.cache.data = UInt8[]
        for ibuffer = 1:nbuffers
            cbuffer = read!(io_cdata, Vector{UInt8}(undef, buffersize[ibuffer]))
            io.cache.data = [io.cache.data;decompress(UInt8, cbuffer)]
        end
    end
    mb = length(io.cache.data) / 1_000_000
    mb_compressed = length(cdata) / 1_000_000
    mbps_read = mb_compressed/t_read
    mbps_decompress = mb / t_decompress
    mbps = mb / (t_read + t_decompress)
    @debug "...data read and decompressed (effective: $mbps MB/s; decompression: $mbps_decompress MB/s; read: $mbps_read MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end

function compressed_offsets(io::CSeis{T,N,LeftJustifyCompressor}, extentindex) where {T,N}
    fmap = foldmap(io, extentindex)
    nframes = length(fmap)
    _headerlength = headerlength(io)

    offset = 1 + nframes*sizeof(Int)

    hdr_offsets = Vector{Int}(undef, nframes)
    for iframe = 1:nframes
        hdr_offsets[iframe] = iframe == 1 ? (offset) : (hdr_offsets[iframe-1] + fmap[iframe-1]*_headerlength)
    end

    offset = hdr_offsets[end] + fmap[end]*_headerlength

    trc_offsets = Vector{Int}(undef, nframes)
    for iframe = 1:nframes
        trc_offsets[iframe] = iframe == 1 ? offset : (trc_offsets[iframe-1] + fmap[iframe-1]*size(io, 1)*sizeof(io.traceformat))
    end

    trc_offsets,hdr_offsets,fmap,nframes,_headerlength
end

function cache_from_file!(io::CSeis{T,N,LeftJustifyCompressor}, extentindex, regularize) where {T,N}
    if regularize
        @debug "reading and decompressing extent $extentindex..."
        t_read = @elapsed begin
            io.cache.data = Vector{UInt8}(undef, cachesize(io, extentindex))
            cdata = unsafe_wrap(Array, pointer(io.cache.data), (filesize(io.extents[extentindex].container, io.extents[extentindex].name),); own=false)
            read!(io.extents[extentindex].container, io.extents[extentindex].name, cdata)
            io.cache.data[1:length(cdata)] .= cdata
        end

        t_decompress = @elapsed begin
            compressed_trc_offsets,compressed_hdr_offsets,fmap,nframes,_headerlength = compressed_offsets(io, extentindex)

            for iframe = nframes:-1:1
                copyto!(io.cache.data, trcsoffset(io, true, extentindex, io.extents[extentindex].frameindices[iframe]), io.cache.data, compressed_trc_offsets[iframe], fmap[iframe]*size(io,1)*sizeof(io.traceformat))
            end

            for iframe = nframes:-1:1
                copyto!(io.cache.data, hdrsoffset(io, true, extentindex, io.extents[extentindex].frameindices[iframe]), io.cache.data, compressed_hdr_offsets[iframe], fmap[iframe]*_headerlength)
            end

            # regularize the indivdual frames
            for (jframe,iframe) = enumerate(io.extents[extentindex].frameindices)
                regularize!(io, fmap[jframe], getframetrcs(io, true, extentindex, iframe), getframehdrs(io, true, extentindex, iframe))
            end
        end
        mb = length(io.cache.data) / 1_000_000
        mb_compressed = length(cdata) / 1_000_000
        mbps_read = mb_compressed/t_read
        mbps_decompress = mb / t_decompress
        mbps = mb / (t_read + t_decompress)
        @debug "...data read and decompressed (effective: $mbps MB/s; decompression: $mbps_decompress MB/s; read: $mbps_read MB/s -- $mb MB; $mb_compressed compressed MB)"
    else
        notacompressor_cache_from_file!(io, extentindex)
    end
    nothing
end

function cache!(io::CSeis, extentindex::Integer, regularize=true, force=false)
    cachetype = regularize ? CACHE_ALL : CACHE_ALL_LEFT_JUSTIFY
    if extentindex == io.cache.extentindex && io.cache.type == cachetype && !force
        return extentindex
    end

    if io.mode == "w" || io.mode == "r+"
        flush(io)
    end

    if isfile(io.extents[extentindex].container, io.extents[extentindex].name)
        cache_from_file!(io, extentindex, regularize)
    else
        @debug "creating cache..."
        io.cache.data = zeros(UInt8, cachesize(io, extentindex))
        @debug "...done creating cache."
    end
    io.cache.extentindex = extentindex
    io.cache.type = cachetype

    extentindex
end

function cache!(io::CSeis, idx::CartesianIndex, regularize=true, force=false)
    extentindex = extentindex_from_frameindex(io, idx)
    cache!(io, extentindex, regularize, force)
end

cache(io::CSeis) = io.cache

partialcache_error() = error("partial caching is only allowed in read-only \"r\" mode")

function cache_foldmap!(io::Union{CSeis{T,N,NotACompressor}, CSeis{T,N,LeftJustifyCompressor}}, extentindex::Integer, force=false) where {T,N}
    io.mode == "r" || partialcache_error()

    if extentindex == io.cache.extentindex && io.cache.type ∈ (CACHE_FOLDMAP,CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY) && !force
        return extentindex
    end

    if isfile(io.extents[extentindex].container, io.extents[extentindex].name)
        @debug "reading foldmap for extent $extentindex from block-storage..."
        t = @elapsed begin
            io.cache.data = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, cachesize_foldmap(io, extentindex)))
        end
        mb = length(io.cache.data)/1_000_000
        mbps = mb/t
        @debug "...foldmap read ($mbps MB/s -- $mb MB)"
    else
        @debug "creating cache..."
        io.cache.data = zeros(UInt8, cachesize_foldmap(io, extentindex))
        @debug "...done creating cache."
    end
    io.cache.extentindex = extentindex
    io.cache.type = CACHE_FOLDMAP

    extentindex
end

function cache_foldmap!(io::CSeis{T,N,BloscCompressor}, extentindex::Integer, force=false) where {T,N}
    if extentindex == io.cache.extentindex && io.cache.type ∈ (CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY) && !force
        return extentindex
    end
    # TODO - we can't to a partial cache here because of how the compression works.
    cache!(io, extentindex, true, force)
    extentindex
end

function cache_foldmap!(io::CSeis, idx::CartesianIndex, force=false)
    extentindex = extentindex_from_frameindex(io, idx)
    cache_foldmap!(io, extentindex, force)
end

function Base.flush(io::CSeis{T,N,NotACompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end
    @debug "writing extent $(io.cache.extentindex) to block-storage..."
    t = @elapsed write(io.extents[io.cache.extentindex].container, io.extents[io.cache.extentindex].name, io.cache.data)
    mb = length(io.cache.data)/1_000_000
    mbps = mb/t
    @debug "...data wrote ($mbps MB/s -- $mb MB)"
    nothing
end

function Base.flush(io::CSeis{T,N,BloscCompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end

    function flush_blosc_buffer_byterange(ibuffer, buffer_size, buffer_remainder)
        isremainder = ibuffer <= buffer_remainder
        firstbyte = (ibuffer - 1)*buffer_size + (isremainder ? ibuffer : buffer_remainder + 1)
        lastbyte = firstbyte + (isremainder ? buffer_size : buffer_size - 1)
        firstbyte,lastbyte
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        Blosc.set_num_threads(Sys.CPU_THREADS)
        Blosc.set_compressor(io.cache.compressor.algorithm)

        maxbuffersize = 2_000_000_000 # This is due to a limitation of the Blosc library.
        cachesize = length(io.cache.data)
        nbuffers,remaining_bytes = divrem(cachesize, maxbuffersize)
        remaining_bytes > 0 && (nbuffers += 1)
        buffersize,buffer_remainder = divrem(cachesize, nbuffers)

        cdata = zeros(UInt8, 8*(1 + nbuffers)) # store number of buffers and length of each compressed buffer
        cdata_io = IOBuffer(cdata; read=true, write=true)
        write(cdata_io, nbuffers)
        close(cdata_io)
        for ibuffer = 1:nbuffers
            firstbyte,lastbyte = flush_blosc_buffer_byterange(ibuffer, buffersize, buffer_remainder)
            _data = unsafe_wrap(Array, pointer(io.cache.data)+(firstbyte-1), (lastbyte-firstbyte+1,); own=false)
            _cdata = compress(_data; level=io.cache.compressor.level, shuffle=io.cache.compressor.shuffle)
            cdata = [cdata;_cdata]
            cdata_io = IOBuffer(cdata; read=true, write=true)
            seek(cdata_io, 8*ibuffer)
            write(cdata_io, length(_cdata))
            close(cdata_io)
        end
        close(cdata_io)
    end
    t_write = @elapsed write(io.extents[io.cache.extentindex].container, io.extents[io.cache.extentindex].name, cdata)
    mb = length(io.cache.data)/1_000_000
    mb_compressed = length(cdata)/1_000_000
    mbps_write = mb_compressed/t_write
    mbps_compress = mb/t_compress
    mbps = mb/(t_write+t_compress)
    @debug "...data compressed and wrote (effective: $mbps MB/s ; compression: $mbps_compress MB/s ; write: $mbps_write MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end

function Base.flush(io::CSeis{T,N,LeftJustifyCompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        for iframe in io.extents[io.cache.extentindex].frameindices
            fmap = foldmap(io)
            fld = fmap[iframe - io.extents[io.cache.extentindex].frameindices[1] + 1]
            leftjustify!(io, getframetrcs(io, true, io.cache.extentindex, iframe), getframehdrs(io, true, io.cache.extentindex, iframe), fld)
        end

        compressed_trc_offsets,compressed_hdr_offsets,fmap,nframes,_headerlength = compressed_offsets(io, io.cache.extentindex)

        for (iframe, jframe) in enumerate(io.extents[io.cache.extentindex].frameindices)
            copyto!(io.cache.data, compressed_hdr_offsets[iframe], io.cache.data, hdrsoffset(io, true, io.cache.extentindex, jframe), fmap[iframe]*_headerlength)
        end

        for (iframe, jframe) in enumerate(io.extents[io.cache.extentindex].frameindices)
            copyto!(io.cache.data, compressed_trc_offsets[iframe], io.cache.data, trcsoffset(io, true, io.cache.extentindex, jframe), fmap[iframe]*size(io,1)*sizeof(io.traceformat))
        end
    end
    t_write = @elapsed begin
        cdata = unsafe_wrap(Array, pointer(io.cache.data), (compressed_trc_offsets[end] + fmap[end]*size(io,1)*sizeof(io.traceformat) - 1,); own=false)
        write(io.extents[io.cache.extentindex].container, io.extents[io.cache.extentindex].name, cdata)
    end

    mb = length(io.cache.data)/1_000_000
    mb_compressed = length(cdata)/1_000_000
    mbps_write = mb_compressed/t_write
    mbps_compress = mb/t_compress
    mbps = mb/(t_write+t_compress)
    @debug "...data compressed and wrote (effective: $mbps MB/s ; compression: $mbps_compress MB/s ; write: $mbps_write MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end

"""
    prop(io::CSeis, propertydef)

Given a property definition `propertydef`, return a trace property corresponding
to the CloudSeis data-set.  `propertydef` can either be a string or a `TracePropertyDef`.
"""
TeaSeis.prop(io::CSeis, _property::Symbol) = io.traceproperties[_property]
TeaSeis.prop(io::CSeis, _property::String) = prop(io, Symbol(_property))
TeaSeis.prop(io::CSeis, _property::Symbol, _T::Type{T}) where {T} = io.traceproperties[_property]::TraceProperty{T}
TeaSeis.prop(io::CSeis, _property::String, _T::Type{T}) where {T} = io.traceproperties[Symbol(_property)]::TraceProperty{T}
TeaSeis.prop(io::CSeis, _property::TracePropertyDef{T}) where {T} = prop(io, _property.label, T)

function foldmap(io::CSeis, extentindex)
    nbytes = sizeof(Int)*length(io.extents[extentindex].frameindices)
    reinterpret(Int, view(io.cache.data, 1:nbytes))
end
foldmap(io::CSeis) = foldmap(io, io.cache.extentindex)

unsafe_foldmap(io::CSeis, extentindex) = unsafe_wrap(Array, convert(Ptr{Int}, pointer(io.cache.data)), (length(io.extents[extentindex].frameindices)); own=false)
unsafe_foldmap(io::CSeis) = unsafe_foldmap(io, io.cache.extentindex)

"""
    fold(io::CSeis, i)

Return the fold of a CloudSeis frame corresponding, and where the frame
is described by its index `i`.  The frame index can be either an
`Int`, a `VarArg{Int}`, or a `CartesianIndex`.

# Example
```julia
io = csopen("foo.cs", "w"; axis_lengths=[10,11,12])
fold(io, 5, 6) # fold corresponding to 5th frame and 6th volume
```
"""
function TeaSeis.fold(io::CSeis, idx::CartesianIndex)
    extentindex = io.mode == "r" ? cache_foldmap!(io, idx) : cache!(io, idx)
    fmap = foldmap(io)
    fmap[linearframeidx(io, idx) - io.extents[extentindex].frameindices[1] + 1]
end
TeaSeis.fold(io::CSeis, idx...) = fold(io, CartesianIndex(idx))

function TeaSeis.fold(io::CSeis, hdrs::AbstractArray{UInt8,2})
    fld = 0
    proptype = prop(io, stockprop[:TRC_TYPE])
    for i = 1:size(hdrs,2)
        if get(proptype, hdrs, i) == tracetype[:live]
            fld += 1
        end
    end
    fld
end

function fold!(io::CSeis, fld, idx::CartesianIndex)
    _foldmode_error() = error("unable to set fold in read-only mode")
    io.mode == "r" && _foldmode_error()
    extentindex = cache!(io, idx)
    fmap = foldmap(io)
    fmap[linearframeidx(io, idx) - io.extents[extentindex].frameindices[1] + 1] = fld
    nothing
end
fold!(io::CSeis, fld, idx...) = fold!(io, fld, CartesianIndex(idx))

"""
    regularize!(io::CSeis[, fld], trcs, hdrs)

This is primarily used for compression, but one can also use this as a convenience method
undo the effect of `leftjustify!`, moving traces from the beginning of the frame to their
recorded trace location.
"""
function TeaSeis.regularize!(io::CSeis, fld, trcs::AbstractArray{Float32, 2}, hdrs::AbstractArray{UInt8, 2})
    if fld == io.axis_lengths[2]
        return
    end
    nsamp, nhead, ntrcs = size(trcs,1), size(hdrs,1), size(trcs,2)
    proptrc, proptyp = prop(io, io.axis_propdefs[2]), prop(io, stockprop[:TRC_TYPE])
    trace_mask = zeros(Int, ntrcs)
    for i = fld:-1:1
        ii = lineartraceidx(io, get(proptrc, hdrs, i))
        trace_mask[ii] = 1
        trcs[:,ii] .= trcs[:,i]
        hdrs[:,ii] .= hdrs[:,i]
    end
    for i = 1:ntrcs
        if trace_mask[i] == 0
            set!(proptrc, hdrs, i, logicaltraceidx(io, i))
            set!(proptyp, hdrs, i, tracetype[:dead])
            trcs[:,i] .= 0
        end
    end
end
TeaSeis.regularize!(io::CSeis, trcs::AbstractArray{Float32, 2}, hdrs::AbstractArray{UInt8, 2}) = regularize!(io, fold(io, hdrs), trcs, hdrs)

"""
    leftjustify!(io::CSeis, trcs, hdrs)

This is primarily used for compression, but one can also use this as a convenience method
to move all live traces to the beginning of its frame.
"""
function TeaSeis.leftjustify!(io::CSeis, trcs::AbstractArray{Float32, 2}, hdrs::AbstractArray{UInt8, 2}, fld::Integer)
    if fld == io.axis_lengths[2]
        return
    end
    proptyp = prop(io, "TRC_TYPE", Int32)
    j,jₒ,ntrcs,nsamp,nhead = 1,2,size(trcs,2),size(trcs,1),size(hdrs,1)
    for i = 1:fld
        if get(proptyp, hdrs, i) != tracetype[:live]
            jₒ = max(jₒ, i + 1)
            for j = jₒ:ntrcs
                if get(proptyp, hdrs, j) == tracetype[:live]
                    @inbounds for k = 1:nsamp
                        tmp = trcs[k,i]
                        trcs[k,i] = trcs[k,j]
                        trcs[k,j] = tmp
                    end

                    @inbounds for k = 1:nhead
                        tmp = hdrs[k,i]
                        hdrs[k,i] = hdrs[k,j]
                        hdrs[k,j] = tmp
                    end
                    jₒ = j + 1
                    break
                end
            end
        end
    end
end
TeaSeis.leftjustify!(io::CSeis, trcs::AbstractArray{Float32, 2}, hdrs::AbstractArray{UInt8, 2}) = leftjustify!(io, trcs, hdrs, fold(io, hdrs))

"""
    leftjustify!(io::CSeis, hdrs)

Convenience method to move all live traces to the beginning of its frame.
"""
function TeaSeis.leftjustify!(io::CSeis, hdrs::AbstractArray{UInt8, 2}, fld::Integer)
    if fld == io.axis_lengths[2]
        return
    end
    proptyp = prop(io, "TRC_TYPE", Int32)
    j,jₒ,ntrcs,nhead = 1,2,size(hdrs, 2),size(hdrs,1)
    for i = 1:fld
        if get(proptyp, hdrs, i) != tracetype[:live]
            jₒ = max(jₒ, i + 1)
            for j = jₒ:ntrcs
                if get(proptyp, hdrs, j) == tracetype[:live]
                    @inbounds for k = 1:nhead
                        tmp = hdrs[k,i]
                        hdrs[k,i] = hdrs[k,j]
                        hdrs[k,j] = tmp
                    end
                    jₒ = j + 1
                    break
                end
            end
        end
    end
end
TeaSeis.leftjustify!(io::CSeis, hdrs::AbstractArray{UInt8, 2}) = leftjustify!(io, hdrs, fold(io, hdrs))

"""
    get(prop::TraceProperty, hdr::Vector)

Returns the vale of a trace header property for the trace header `hdr`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
hdrs = readframehdrs(io,1)
hdr = @view hdrs[:,1]
get(prop(io, "TRACE"), hdr)
close(io)
```
"""
function Base.get(prop::TraceProperty{T}, hdr::AbstractVector{UInt8}) where {T<:Number}
    iohdr = IOBuffer(hdr, read=true)
    seek(iohdr, prop.byteoffset)
    read(iohdr, T)
end

function Base.get(prop::TraceProperty{T}, hdr::AbstractVector{UInt8}) where {T<:Vector}
    iohdr = IOBuffer(hdr, read=true)
    seek(iohdr, prop.byteoffset)
    read!(iohdr, T(undef, prop.def.elementcount))
end

"""
    get(prop::TraceProperty, hdrs::Matrix, i)

Returns the vale of a trace header property for the ith column in the trace headers `hdrs`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
hdrs = readframehdrs(io, 1)
get(prop(io, "TRACE"), hdrs, 1)
close(io)
```
"""
Base.get(prop::TraceProperty, hdrs::AbstractArray{UInt8,2}, i::Integer) = get(prop, @view(hdrs[:,i]))

"""
    set!(prop::TraceProperty, hdr::Vector, value)

Set the vale of a trace header property for the `hdr`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
hdrs = readframehdrs(io, 1)
hdr = @view hdrs[:,1]
set!(prop(io, "TRACE"), hdr, 1)
close(io)
```
"""
function TeaSeis.set!(prop::TraceProperty{T}, hdr::AbstractArray{UInt8,1}, value) where {T}
    iohdr = IOBuffer(hdr, read=true, write=true)
    seek(iohdr, prop.byteoffset)
    write(iohdr, T(value))
    nothing
end

"""
    set!(prop::TraceProperty, hdrs::Matrix, i, value)

Set the vale of the header property for the ith column of `hdrs`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
hdrs = readframehdrs(io, 1)
set!(prop(io, "TRACE"), hdr, 1, 2)
close(io)
```
"""
TeaSeis.set!(prop::TraceProperty, hdrs::AbstractArray{UInt8,2}, i::Integer, value) = set!(prop, @view(hdrs[:,i]), value)

Base.LinearIndices(io::CSeis) = LinearIndices(size(io)[3:end])

function linearframeidx(io, idx, idim)
    @boundscheck begin
        if rem(idx[idim] - io.axis_lstarts[2+idim], io.axis_lincs[2+idim]) > 0
            error("$idx is out-of-bounds.")
        end
    end
    div(idx[idim] - io.axis_lstarts[2+idim], io.axis_lincs[2+idim]) + 1
end

linearframeidx(io, idx::NTuple{N,Int}) where {N} = LinearIndices(io)[CartesianIndex(ntuple(idim->linearframeidx(io, idx, idim), N))]
linearframeidx(io, idx::CartesianIndex) = linearframeidx(io, idx.I)
linearframeidx(io, idx::Int...) = linearframeidx(io, idx)

logicalframeidx(io, idx::NTuple{N,Int}) where {N} = CartesianIndex(ntuple(idim->io.axis_lstarts[2+idim] + io.axis_lincs[2+idim]*(idx[idim] - 1), N))
logicalframeidx(io, idx::CartesianIndex) = logicalframeidx(io, idx.I)
logicalframeidx(io, idx::Int...) = logicalframeidx(io, idx)

function lineartraceidx(io, idx)
    @boundscheck begin
        if rem(idx - io.axis_lstarts[2], io.axis_lincs[2]) > 0
            error("$idx is out-of-bounds.")
        end
    end
    div(idx - io.axis_lstarts[2], io.axis_lincs[2]) + 1
end

logicaltraceidx(io, idx) = io.axis_lstarts[2] + io.axis_lincs[2]*(idx-1)

"""
    I = LogicalIndices(io::CSeis)

Returns a construct similar to `CartesianIndices` that allows conversion from linear
indices to cartesian indices that are offest by the logical starts and deltas of the
CloudSeis data context.  In addition, `LogicalIndices` implements iteration for looping
over all frames in a data-set.  For example,

```julia
io = csopen(AzContainer("foo"; storageaccount="bar))
idx = LogicalIndices(io)[2] # get the index corresponding to the second frame in the data-set.
for idx in LogicalIndices(io)
    @show idx
    trcs, hdrs = getframe(io, idx)
end
```
"""
struct LogicalIndices{T,N}
    io::CSeis{T,N}
end

Base.getindex(c::LogicalIndices, i) = logicalframeidx(c.io, CartesianIndices(c.io.axis_lengths[3:end])[i])
Base.length(c::LogicalIndices) = prod(size(c.io)[3:end])

Base.iterate(c::LogicalIndices, state=1) = state > length(c) ? nothing : (c[state], state+1)

function extentindex_from_frameindex(io, idx::CartesianIndex)
    frameidx = linearframeidx(io, idx)

    # TODO - optimize me
    for (i,extent) in enumerate(io.extents)
        if frameidx ∈ extent.frameindices
            return i
        end
    end
    error("can't find extent")
end

function trcsoffset(io::CSeis, regularize::Bool, extentindex, frameidx)
    nsamples = size(io,1)
    nframes = length(io.extents[extentindex].frameindices)
    frstframe = io.extents[extentindex].frameindices[1]

    local ntraces_times_nframes, ntraces_times_frameidx
    if regularize
        ntraces = size(io,2)
        ntraces_times_nframes = ntraces*nframes
        ntraces_times_frameidx = ntraces*(frameidx-frstframe)
    else
        ntraces_times_nframes = mapreduce(frameindex->fold(io,frameindex), +, io.extents[extentindex].frameindices)
        ntraces_times_frameidx = frstframe == frameidx ? 0 : mapreduce(frameidx->fold(io,frameidx), +, frstframe:(frameidx-1))
    end

    sizeof(Int)*nframes + io.hdrlength*ntraces_times_nframes + sizeof(io.traceformat)*nsamples*ntraces_times_frameidx + 1
end

function hdrsoffset(io::CSeis, regularize::Bool, extentindex, frameidx)
    ntraces = size(io,2)
    nframes = length(io.extents[extentindex].frameindices)
    frstframe = io.extents[extentindex].frameindices[1]

    local ntraces_times_frameidx
    if regularize
        ntraces_times_frameidx = ntraces*(frameidx-frstframe)
    else
        ntraces_times_frameidx = frstframe == frameidx ? 0 : mapreduce(_frameidx->fold(io,_frameidx), +, frstframe:(frameidx-1))
    end

    offset = sizeof(Int)*nframes + io.hdrlength*ntraces_times_frameidx + 1
    offset
end

function unsafe_gethdrs(io::CSeis, extentindex)
    regularize = true # TODO

    frameindices = io.extents[extentindex].frameindices

    byteoffset = hdrsoffset(io, regularize, extentindex, frameindices[1]) - 1

    ntraces = size(io, 2)
    nframes = length(frameindices)

    rem(io.hdrlength,4) == 0 || error("header length must be a scalar multiple of 4")

    unsafe_wrap(Array, convert(Ptr{Int32}, pointer(io.cache.data) + byteoffset), (div(io.hdrlength,4),ntraces,nframes); own=false)
end
unsafe_gethdrs(io::CSeis) = unsafe_gethdrs(io, io.cache.extentindex)

function unsafe_gettrcs(io::CSeis{T}, extentindex) where {T}
    regularize = true # TODO

    frameindices = io.extents[extentindex].frameindices

    byteoffset = trcsoffset(io, regularize, extentindex, frameindices[1]) - 1

    ntraces = size(io, 2)
    nframes = length(frameindices)

    unsafe_wrap(Array, convert(Ptr{T}, pointer(io.cache.data) + byteoffset), (size(io,1), ntraces, nframes), own=false)
end
unsafe_gettrcs(io::CSeis) = unsafe_gettrcs(io, io.cache.extentindex)

function getframetrcs(io::CSeis, regularize::Bool, extentindex, idx::Int)
    data = io.cache.data

    ntraces = regularize ? size(io,2) : fold(io,idx)

    if ntraces == 0
        return io.traceformat[]
    end

    frstbyte = trcsoffset(io, regularize, extentindex, idx)
    lastbyte = frstbyte + ntraces*size(io,1)*sizeof(io.traceformat) - 1
    reshape(reinterpret(io.traceformat, view(data,frstbyte:lastbyte)), :, ntraces)
end
getframetrcs(io::CSeis, regularize::Bool, idx::CartesianIndex) = getframetrcs(io, regularize, cache!(io, idx, regularize), linearframeidx(io, idx))
getframetrcs(io::CSeis, regularize::Bool, idx...) = getframetrcs(io, regularize, CartesianIndex(idx))

function getframehdrs(io::CSeis, regularize::Bool, extentindex, idx::Int)
    data = io.cache.data

    ntraces = regularize ? size(io,2) : fold(io,idx)

    if ntraces == 0
        return UInt8[]
    end

    frstbyte = hdrsoffset(io, regularize, extentindex, idx)
    lastbyte = frstbyte + ntraces*io.hdrlength - 1
    reshape(view(data,frstbyte:lastbyte), :, ntraces)
end
getframehdrs(io::CSeis, regularize::Bool, idx::CartesianIndex) = getframehdrs(io, regularize, cache!(io, idx, regularize), linearframeidx(io, idx))
getframehdrs(io::CSeis, regularize::Bool, idx...) = getframehdrs(io, regularize, CartesianIndex(idx))

getframe(io::CSeis, regularize::Bool, idx) = getframetrcs(io, regularize, idx), getframehdrs(io, regularize, idx)

"""
    trcs = allocframetrcs(io::CSeis)

Allocate and return memory to store traces for a single frame.
"""
TeaSeis.allocframetrcs(io::CSeis) = zeros(Float32, size(io,1), size(io,2))

"""
    hdrs = allocframehdrs(io::CSeis)

Allocate and return memory to store trace headers for a single frame.
"""
TeaSeis.allocframehdrs(io::CSeis) = zeros(UInt8, headerlength(io), size(io,2))

"""
    trcs,hdrs = allocframe(io::CSeis)

Allocate and return memory to store traces and trace headers for a single frame.
"""
TeaSeis.allocframe(io::CSeis) = allocframetrcs(io), allocframehdrs(io)

function regularization_check(io, regularize)
    if !regularize && (!isa(io.cache.compressor, LeftJustifyCompressor) || io.mode != "r")
        error("regularize=false is only valid for compression method 'LeftJustifyCompressor' and mode read only 'r'.")
    end
end

"""
    readframetrcs!(io::CSeis, trcs, idx...[; regularize=true])

Read traces from `io` into `trcs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
function TeaSeis.readframetrcs!(io::CSeis, trcs::AbstractArray, idx::CartesianIndex; regularize=true)
    @boundscheck regularization_check(io, regularize)
    _trcs = getframetrcs(io, regularize, idx)
    copyto!(trcs, 1, _trcs, 1, length(_trcs))
end

"""
    trcs = readframetrcs(io::CSeis, idx...; regularize=true)

Read traces from `io` for the frame `idx...`.  `idx..` can either be
integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
TeaSeis.readframetrcs(io::CSeis, idx::CartesianIndex; regularize=true) = readframetrcs!(io, allocframetrcs(io), idx; regularize)

TeaSeis.readframetrcs!(io::CSeis, trcs::AbstractArray, idx...; regularize=true) = readframetrcs!(io, trcs, CartesianIndex(idx); regularize)
TeaSeis.readframetrcs(io::CSeis, idx...; regularize=true) = readframetrcs(io, CartesianIndex(idx); regularize)

"""
    readframehdrs!(io::CSeis, trcs, idx...; regularize=true)

Read headers from `io` into `hdrs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
function TeaSeis.readframehdrs!(io::CSeis, hdrs::AbstractArray{UInt8,2}, idx::CartesianIndex; regularize=true)
    @boundscheck regularization_check(io, regularize)
    _hdrs = getframehdrs(io, regularize, idx)
    copyto!(hdrs, 1, _hdrs, 1, length(_hdrs))
    fld = fold(io, idx)
    if !regularize || fld == 0
        prop_trctype = prop(io, stockprop[:TRC_TYPE])
        for itrace = (fld+1):size(io, 2)
            set!(prop_trctype, hdrs, itrace, tracetype[:dead])
        end
    end
    hdrs
end

"""
    hdrs = readframehdrs(io::CSeis, idx...; regularize=true)

Read headers from `io` for the frame `idx...`.  `idx..` can either be
integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
TeaSeis.readframehdrs(io::CSeis, idx::CartesianIndex; regularize=true) = readframehdrs!(io, allocframehdrs(io), idx; regularize)

TeaSeis.readframehdrs!(io::CSeis, hdrs::AbstractArray{UInt8,2}, idx...; regularize=true) = readframehdrs!(io, hdrs, CartesianIndex(idx); regularize)
TeaSeis.readframehdrs(io::CSeis, idx...; regularize=true) = readframehdrs(io, CartesianIndex(idx); regularize)

"""
    readframe!(io::CSeis, trcs, hdrs, idx...; regularize=true)

Read traces and headers from `io` into `trcs::Matrix`, and `hdrs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
TeaSeis.readframe!(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2}, idx::CartesianIndex; regularize=true) = readframetrcs!(io, trcs, idx; regularize), readframehdrs!(io, hdrs, idx; regularize)

"""
    trcs,hdrs = readframe(io::CSeis, idx...; regularize=true)

Read traces and headers from `io` for the frame `idx...`.  `idx...` can
either be integer(s) or a `CartesianIndex`.

The `regularize` named argument is only applicable when the compression
method is `LeftJustifyCompressor`.  If set to true, then traces are 
regularized to their correct context locations.  Otherwise, they remain
left justified.  Note that one can subsequently use the `regularize!` method.
"""
TeaSeis.readframe(io::CSeis, idx::CartesianIndex; regularize=true) = readframe!(io, allocframetrcs(io), allocframehdrs(io), idx; regularize)
TeaSeis.readframe!(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2}, idx...; regularize=true) = readframe!(io, trcs, hdrs, CartesianIndex(idx); regularize)
TeaSeis.readframe(io::CSeis, idx...; regularize=true) = readframe(io, CartesianIndex(idx); regularize)

function get_first_live_trace(io, hdrs)
    prop_trctype = prop(io, stockprop[:TRC_TYPE])
    for itrace = 1:size(io,2)
        if get(prop_trctype, hdrs, itrace) == tracetype[:live]
            return itrace
        end
    end
    return size(io,2)+1
end

"""
    writeframe(io::CSeis, trcs, hdrs)

Write a frame to `io`.  The location of the frame is determined
by the axis headers set in `hdrs`.
"""
function TeaSeis.writeframe(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2})
    j = get_first_live_trace(io, hdrs)
    if j > size(io, 2)
        return 0
    end
    idx = CartesianIndex(ntuple(i->get(prop(io, io.axis_propdefs[i+2]), hdrs, j), ndims(io)-2))
    cache!(io, idx)
    fld = fold(io, hdrs)
    fold!(io, fld, idx)
    _trcs = getframetrcs(io, true, idx)
    _hdrs = getframehdrs(io, true, idx)
    copyto!(_trcs, 1, trcs, 1, size(io,2)*size(io,1))
    copyto!(_hdrs, 1, hdrs, 1, size(io,2)*io.hdrlength)
    fld
end

"""
    writeframe(io::CSeis, trcs, idx...)

Write a frame to `io`.  The location of the frame is determined from
`idx...` which can either be integer(s) or a `CartesianIndex`.  A minimal
set of headers are created from `idx...` and are also written to `io`.
""" 
function TeaSeis.writeframe(io::CSeis, trcs::AbstractArray, idx::CartesianIndex)
    cache!(io, idx)
    _trcs,hdrs = getframe(io, true, idx)
    prop_trctype = prop(io, stockprop[:TRC_TYPE])
    for i = 1:size(io,2)
        set!(prop_trctype, hdrs, i, tracetype[:live])
        set!(prop(io, io.axis_propdefs[2]), hdrs, i, i)
        for idim = 3:ndims(io)
            set!(prop(io, io.axis_propdefs[idim]), hdrs, i, idx[idim-2])
        end
    end
    copyto!(_trcs, 1, trcs, 1, size(io,1)*size(io,2))
    fold!(io, size(io,2), idx)
    size(io,2)
end
TeaSeis.writeframe(io::CSeis, trcs::AbstractArray, idx...) = writeframe(io, trcs, CartesianIndex(idx))

function parserngs(io::CSeis, smprng::Union{Int,AbstractRange{Int},Colon}, trcrng::Union{Int,AbstractRange{Int},Colon}, rng::Vararg{Union{Int,AbstractRange{Int},Colon},N}) where N
    smprng = parserng(io, smprng, 1)
    trcrng = parserng(io, trcrng, 2)
    rng = ntuple(i->parserng(io, rng[i], 2+i), N)
    nrng = ntuple(i->length(rng[i]), N)::NTuple{N,Int}
    smprng::StepRange{Int,Int}, trcrng::StepRange{Int,Int}, rng::NTuple{N,StepRange{Int,Int}}, nrng::NTuple{N,Int}
end
parserng(io::CSeis,rng::Int,i) = StepRange(rng,1,rng)
parserng(io::CSeis,rng::StepRange{Int},i) = rng
parserng(io::CSeis,rng::AbstractRange{Int},i) = StepRange(rng)
parserng(io::CSeis,rng::Colon,i) = (1:1:size(io,i))::StepRange{Int,Int}

parseindex(rng::Tuple{AbstractRange}, idx_n::CartesianIndex{1}) = CartesianIndex(rng[1][idx_n[1]])
parseindex(rng::Tuple{AbstractRange,AbstractRange}, idx_n::CartesianIndex{2}) = CartesianIndex(rng[1][idx_n[1]],rng[2][idx_n[2]])
parseindex(rng::Tuple{AbstractRange,AbstractRange,AbstractRange}, idx_n::CartesianIndex{3}) = CartesianIndex(rng[1][idx_n[1]],rng[2][idx_n[2]],rng[3][idx_n[3]])
parseindex(rng::NTuple{N,AbstractRange}, idx_n::CartesianIndex) where {N} = CartesianIndex(ntuple(i->rng[i][idx_n[i]], length(rng)))

function readtrcs_impl!(io::CSeis, trcs::AbstractArray, smprng::AbstractRange{Int}, trcrng::AbstractRange{Int}, rng::Vararg{Union{Colon,Int,AbstractRange{Int}},N}) where {N}
    n = ntuple(i->length(rng[i]), N)::NTuple{N,Int}
    for idx_n in CartesianIndices(n)
        idx = logicalframeidx(io, parseindex(rng, idx_n))
        frmtrcs = getframetrcs(io, true, idx)
        for (itrc,trc) in enumerate(trcrng), (ismp,smp) in enumerate(smprng)
            trcs[ismp,itrc,idx_n] = frmtrcs[smp,trc]
        end
    end
    trcs
end

"""
    readtrcs!(io::CSeis, trcs, rng...)

Read traces from `io` into `trcs::AbstractArray`, and where the size and dimension
of `trcs` must correspond to the range specification in `rng...`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
size(io) # (100,101,102)
trcs = readtrcs!(io, Array{Float32,3}(undef,50,101,1), 1:50, :, 1)
close(io)
```
"""
function TeaSeis.readtrcs!(io::CSeis, trcs::AbstractArray, rng::Vararg{Union{Int,AbstractRange{Int},Colon},N}) where {N}
    smprng, trcrng, _rng, nrng = parserngs(io, rng...)
    @assert size(trcs)[1:2] == (length(smprng), length(trcrng)) && size(trcs)[3:end] == nrng
    readtrcs_impl!(io, trcs, smprng, trcrng, _rng...)
end

"""
    readtrcs(io::CSeis, trcs, rng...)

Read traces from `io` for the specified `rng...`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
size(io) # (100,101,102)
trcs = readtrcs(io, 1:50, :, 1)
close(io)
```
"""
function TeaSeis.readtrcs(io::CSeis, rng::Vararg{Union{Int,AbstractRange{Int},Colon},N}) where {N}
    smprng, trcrng, _rng, nrng = parserngs(io, rng...)
    trcs = zeros(Float32, length(smprng), length(trcrng), nrng...)
    readtrcs_impl!(io, trcs, smprng, trcrng, _rng...)
end

function readhdrs_impl!(io::CSeis, hdrs::AbstractArray, trcrng::AbstractRange{Int}, rng::Vararg{Union{Colon,Int,AbstractRange{Int}},N}) where {N}
    n = ntuple(i->length(rng[i]), N)::NTuple{N,Int}
    for idx_n in CartesianIndices(n)
        idx = logicalframeidx(io, parseindex(rng, idx_n))
        frmhdrs = getframehdrs(io, true, idx)
        for (itrc,trc) in enumerate(trcrng), ismp in 1:headerlength(io)
            hdrs[ismp,itrc,idx_n] = frmhdrs[ismp,trc]
        end
    end
    hdrs
end

"""
    readhdrs!(io::CSeis, hdrs, rng...)

Read trace headers from `io` into `hdrs::AbstractArray`, and where the size and dimension
of `hdrs` must correspond to the range specification in `rng...`.  Unlike `readtrcs!`, we
always read the entirity of each trace header.  Therefore, the first range specifier must
be `:`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
size(io) # (100,101,102)
hdrs = readhdrs!(io, Array{Float32,3}(undef,headerlength(io),101,1), :, :, 1)
close(io)
```
"""
function TeaSeis.readhdrs!(io::CSeis, hdrs::AbstractArray, rng::Vararg{Union{Int,AbstractRange{Int},Colon},N}) where {N}
    smprng, trcrng, _rng, nrng = parserngs(io, rng...)
    @assert size(hdrs)[1:2] == (headerlength(io), length(trcrng)) && size(hdrs)[3:end] == nrng
    readhdrs_impl!(io, hdrs, smprng, trcrng, _rng...)
end

"""
    readhdrs(io::CSeis, rng...)

Read trace headers from `io`.  Unlike `readtrcs!`, we always read the entirity of
each trace header.  Therefore, the first range specifier must be `:`.

# Example
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container)
size(io) # (100,101,102)
hdrs = readhdrs(io, :, :, 1)
close(io)
```
"""
function TeaSeis.readhdrs(io::CSeis, rng::Vararg{Union{Int,AbstractRange{Int},Colon},N}) where {N}
    smprng, trcrng, _rng, nrng = parserngs(io, rng...) # smprng is not used
    hdrs = zeros(UInt8, headerlength(io), length(trcrng), nrng...)
    readhdrs_impl!(io, hdrs, trcrng, _rng...)
end

"""
    write(io::CSeis, trcs, rng...)

Write traces to `io`, and where the size and dimension of `trcs` must match
the size and dimensions of `rng...`.

# Example 1
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container, "w"; axis_length=[100,101,102])
write(io, rand(Float32,100,101,102), :, :, :)
close(io)
```

# Example 2
```julia
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container, "w"; axis_lengths=[100,101,102])
write(io, rand(Float32,100,50,1), :, 1:50, 1)
close(io)
```
"""
function Base.write(io::CSeis, trcs::AbstractArray, rng::Vararg{Union{Colon,Int,AbstractRange{Int}},N}) where {N}
    smprng, trcrng, _rng, nrng = parserngs(io, rng...)
    frmtrcs = allocframetrcs(io)
    write_helper(io, trcs, frmtrcs, smprng, trcrng, nrng, _rng)
end

function write_helper(io::CSeis, trcs, frmtrcs, smprng, trcrng, nrng, _rng::NTuple{N,StepRange{Int,Int}}) where {N}
    n = ntuple(i->length(_rng[i]), N)::NTuple{N,Int}

    frmhdrs = allocframehdrs(io)
    prop_trctype = prop(io, stockprop[:TRC_TYPE])
    map(itrace->set!(prop_trctype, frmhdrs, itrace, tracetype[:live]), 1:size(io,2))

    props = map(idim->prop(io,io.axis_propdefs[idim]), 1:ndims(io))
    for idx_n in CartesianIndices(n)
        idx = logicalframeidx(io, parseindex(_rng, idx_n))

        frmtrcs = getframetrcs(io, true, idx)
        for (itrc,trc) in enumerate(trcrng), (ismp,smp) in enumerate(smprng)
            frmtrcs[smp,trc] = trcs[ismp,itrc,idx_n]
        end

        for itrace = 1:io.axis_lengths[2]
            set!(props[2], frmhdrs, itrace, logicaltraceidx(io, itrace))
            for idim = 3:ndims(io)
                set!(props[idim], frmhdrs, itrace, idx[idim-2])
            end
        end
        writeframe(io, frmtrcs, frmhdrs)
    end
end

"""
    reduce(io[; optional keyword arguments...])

Consolodate data-set extents into fewer extents.

# Optional keyword arguments
* `mbytes_per_extent=1000` set the size of each output extent via either the size of each extent
* `frames_per_extent=0` set the size of each output extent via the number of frames allowed in each extent.
"""
function Base.reduce(io::CSeis; mbytes_per_extent=1000, frames_per_extent=0)
    local nextents
    nframes = prod(size(io)[3:end])
    if frames_per_extent > 0
        nextents,r = divrem(prod(io.axis_lengths[3:end]), frames_per_extent)
        if r > 0
            nextents += 1
        end
    else # mbytes_per_extent
        bytes_total = prod(size(io)[2:end])*(size(io,1)*sizeof(io.traceformat) + headerlength(io)) + sizeof(Int)*nframes
        nextents = max(1, div(bytes_total, mbytes_per_extent*1000*1000))
    end
    nominal_frames_per_extent, remaining_frames = divrem(nframes, nextents)

    r = randstring(4)
    extents = make_extents(io.containers, nextents, "extents-reduced-$r", nominal_frames_per_extent, remaining_frames)

    _io = copy(io, "w", extents)

    function _reduce(tsk, io, _io)
        t,h = allocframe(io)
        I = CartesianIndices(size(io)[3:end])
        for i in _io.extents[tsk].frameindices
            readframe!(io, t, h, I[i])
            writeframe(_io, t, h)
        end
        flush(_io)
        nothing
    end
    pmap(tsk->_reduce(tsk, io, _io), 1:nextents)

    _rm(tsk, io) = rm(io.extents[tsk].container, io.extents[tsk].name)
    pmap(tsk->_rm(tsk, io), 1:length(io.extents))

    description = JSON.parse(read(io.containers[1], "description.json", String))
    description["extents"] = [Dict(extents[i]) for i=1:nextents]
    write(_io.containers[1], "description.json", JSON.json(description, 1))

    _io
end

"""
    empty!(io::CSeis)

Empty (i.e. delete the extents) from a CloudSeis data-set.
"""
function Base.empty!(io::CSeis)
    for container in io.containers
        objects = filter(object->startswith(object, "extent"), readdir(container))
        for object in objects
            rm(container, object)
        end
    end
end

function description_axis_lengths!(io, description, axis_lengths)
    old_axis_lengths = size(io)
    if length(old_axis_lengths) != length(axis_lengths)
        error("mismatch between dataset dimension and 'axis_lengths' parameter")
    end
    for idim = 1:ndims(io)-1
        if old_axis_lengths[idim] != axis_lengths[idim]
            error("unable to mutate the first $(ndims(io)-1) dimensions")
        end
    end
    if prod(axis_lengths[3:end]) < prod(old_axis_lengths[3:end])
        error("'axis_lengths' must contain more frames (2D slices) than the dataset")
    end

    nominal_frames_per_extent = length(io.extents[1].frameindices)
    delta_frames = prod(axis_lengths[3:end]) - prod(old_axis_lengths[3:end])
    nextents = max(div(delta_frames, nominal_frames_per_extent), 1)
    nominal_frames_per_extent,remaining_frames = divrem(delta_frames, nextents)

    foldername = join(split(io.extents[end].name, '/')[1:end-1], '/')

    all_containers = [io.extents[1].container]
    for iextent in 2:length(io.extents)
        io.extents[iextent].container ∈ all_containers && break
        push!(all_containers, io.extents[iextent].container)
    end

    last_container_index = findfirst(container->container == io.extents[end].container, all_containers)
    last_container_index === nothing && error("unable to find the last container index.")

    lastframe = io.extents[end].frameindices[end]

    last_extent_index = 0
    for extent in io.extents
        _last_extent_index = parse(Int, split(splitpath(extent.name)[2], '-')[end])
        if _last_extent_index > last_extent_index
            last_extent_index = _last_extent_index
        end
    end

    new_extents = make_extents(all_containers, nextents, foldername, nominal_frames_per_extent, remaining_frames; last_container_index, last_extent_index, lastframe)
    description["extents"] = [description["extents"]..., Dict.(new_extents)...]
    description["fileproperties"]["axis_lengths"] = axis_lengths
end

"""
    description!(io::CSeis, kwargs...)

Limited mutation of the properties of an existing CloudSeis dataset.

# Optional key-word arguments
* `axis_lengths = nothing`  Grow the size of the axis lengths (3rd dimension and higher)[1].

# Notes
* [1] `axis_lengths` must be the same length as `size(io)`.  In addition, `prod(axis_lengths[3:end])`
must be greater than `prod(size(io)[3:end])` and `axis_lengths[i]` must equal `size(io,i)` for i∈(1..ndims(io)-1).
"""
function description!(io::CSeis; axis_lengths=nothing)
    io.mode == "r+" || error("mutation is only availabe for data-sets open in 'r+' mode.")

    description = JSON.parse(read(joinpath(io.containers[1], "description.json"), String))

    if axis_lengths !== nothing
        description_axis_lengths!(io, description, axis_lengths)
    end

    write(io.containers[1], "description.json", json(description, 1))
end

function cp_extent(iextent, src_extent, dst_extent, nextents)
    if isfile(src_extent.container, src_extent.name)
        _filesize = filesize(src_extent.container, src_extent.name) / 1000 / 1000
        t = @elapsed cp(src_extent.container, src_extent.name, dst_extent.container, dst_extent.name)
        @info "copied extent $iextent/$nextents in $t seconds ($_filesize MB, $(_filesize / t) MB/s)"
    end
    nothing
end

function cp_extents_batch(ibatch, nbatch, batch_size, iextents, src_extents, dst_extents, nextents)
    i1 = (ibatch-1)*batch_size + 1
    i2 = ibatch == nbatch ? length(iextents) : i1 + batch_size - 1
    asyncmap(iextent->cp_extent(iextents[iextent], src_extents[iextent], dst_extents[iextent], nextents), i1:i2)
    nothing
end

"""
    cp(src::CSeis, dst[, extents=:]; batch_size=32, workers=Distributed.workers)

Copy a CloudSeis data-set to `dst` and where `dst` is either of type `Container` or
of type `Vector{Container}`.  The latter is used for sharding data across multiple
storage accounts.

The option `extents` can be used to copy a sub-set of the CloudSeis data extents (e.g. `1:10`).
The `batch_size` option allows extents to be copied in batches and where the number of extents
associated with each batch is set via `batch_size`.  Within a batch, each extent is copied
via an asynchronous task.

The `cp` method will be executed on the set of machines defined by `workers`.  Note that the
work will be sent to the workers one batch at a time.
"""
function Base.cp(src::CSeis, dst_containers::Vector{<:Container}, extents=Colon(); batch_size=32, workers=Distributed.workers)
    description = JSON.parse(read(src.containers[1], "description.json", String))

    src_extents = [Extent(description["extents"][iextent], src.containers) for iextent in eachindex(description["extents"])]
    dst_extents = similar(src_extents)

    for iextent in eachindex(description["extents"])
        k = rem(iextent - 1, length(dst_containers)) + 1
        dst_extents[iextent] = Extent(src_extents[iextent].name, dst_containers[k], src_extents[iextent].frameindices)
        description["extents"][iextent] = Dict(dst_extents[iextent])
    end

    for dst_container in dst_containers
        mkpath(dst_container)
    end

    write(dst_containers[1], "description.json", json(description, 1))

    local iextents
    if extents == Colon()
        iextents = 1:length(description["extents"])
    else
        iextents = extents
    end

    nextents = length(description["extents"])
    nbatch = clamp(div(length(iextents), batch_size), 1, length(iextents))

    pmap(ibatch->cp_extents_batch(ibatch, nbatch, batch_size, iextents, src_extents, dst_extents, nextents), CachingPool(workers()), 1:nbatch)
    nothing
end
Base.cp(src::CSeis, dst_container::Container, extents=Colon(); kwargs...) = cp(src, [dst_container], extents; kwargs...)

"""
    mv(src::CSeis, dst::Container; batch_size=32, workers=Distributed.workers)

move a CloudSeis data-set to `dst`.  See `cp` for a description of `batch_size`
and `workers`.
"""
function Base.mv(src::CSeis, dst::Container; kwargs...)
    cp(src, dst, project=project, location=location; kwargs...)
    rm(src)
end

"""
    ndims(io::CSeis)

Return the number of dimensions in a CloudSeis data-set.
"""
Base.ndims(io::CSeis) = length(io.axis_lengths)

"""
    size(io::CSeis[,i])

Return the size of a CloudSeis data-set.  If `i` is specified,
then return the size of dimension `i`.
"""
Base.size(io::CSeis) = io.axis_lengths
Base.size(io::CSeis,i) = io.axis_lengths[i]

"""
    length(io::CSeis)

Return the total number of frames in a CloudSeis data-set.
"""
Base.length(io::CSeis) = prod(io.axis_lengths[3:end])

"""
    propdefs(io::CSeis[,i])

Return the property definitions corresponding to the data-context axes.
If `i` is specified, then return the property definition for the ith
data-context axis.
"""
TeaSeis.propdefs(io::CSeis) = io.axis_propdefs
TeaSeis.propdefs(io::CSeis,i) = io.axis_propdefs[i]

"""
    props(io::CSeis[,i])

Return the trace property corresponding to the data-context axes.
If `i` is specified, then return the property for the ith data-context
axis.
"""
TeaSeis.props(io::CSeis,i) = prop(io, io.axis_propdefs[i])
TeaSeis.props(io::CSeis) = ntuple(i->props(io, i), ndims(io))

"""
    pincs(io::CSeis[,i])

Return the physical increment for the data-context axes.
If `i` is specified, then return the physical increment for
the ith data-context axis.
"""
TeaSeis.pincs(io::CSeis) = io.axis_pincs
TeaSeis.pincs(io::CSeis,i) = io.axis_pincs[i]

"""
    pstarts(io::CSeis[,i])

Return the physical start for the data-context axes.
If `i` is specified, then return the physical start
for the ith data-context axis.
"""
TeaSeis.pstarts(io::CSeis,i) = io.axis_pstarts[i]
TeaSeis.pstarts(io::CSeis) = io.axis_pstarts

"""
    lincs(io::CSeis[,i])

Return the logical increment for the data-context axes.
If `i` is specfied, then return the logical increment for
the ith data-context axis
"""
TeaSeis.lincs(io::CSeis,i) = io.axis_lincs[i]
TeaSeis.lincs(io::CSeis) = io.axis_lincs

"""
    lstarts(io::CSeis[,i])

Return the logical start for the data-context axes.
If `i` is specfied, then return the logical start for
the ith data-context axis
"""
TeaSeis.lstarts(io::CSeis,i) = io.axis_lstarts[i]
TeaSeis.lstarts(io::CSeis) = io.axis_lstarts

"""
    units(io::CSeis[,i])

Return the physical unit for the data-context axes.
If `i` is specfied, then return the physical unit
for the ith data-context axis.
"""
TeaSeis.units(io::CSeis,i) = io.axis_units[i]
TeaSeis.units(io::CSeis) = io.axis_units

"""
    domains(io::CSeis[,i])

Return the domain for the data-context axes.
If `i` is specified, then return the domain
for the ith data-context axis.
"""
TeaSeis.domains(io::CSeis,i) = io.axis_domains[i]
TeaSeis.domains(io::CSeis) = io.axis_domains

"""
    geometry(io::CSeis)

Return the geometry (if any) associated with the
CloudSeis data-set.
"""
TeaSeis.geometry(io::CSeis) = io.geometry

"""
    in(propdef::TracePropertyDef, io::CSeis)

return `true` if `propdef` exists in the CloudSeis
data-set, `io`.
"""
function Base.in(propdef::TracePropertyDef, io::CSeis)
    for traceprop in io.traceproperties
        if propdef.label == traceprop.def.label
            return true
        end
    end
    false
end

"""
    dataproperty(io::CSeis, name[, default_value])

Return the CloudSeis data property with `name`
that is in the CloudSeis data-set `io`.  If
`default_value` is provided and a property with
`name` does not exist, then return `default_value`.
"""
function dataproperty(io::CSeis, nm)
    try
        return io.dataproperties[Symbol(nm)]
    catch
        error("data property with label $nm does not exist")
    end
end

dataproperty(io::CSeis, nm, default_value) = hasdataproperty(io::CSeis, nm) ? dataproperty(io, nm) : default_value

"""
    hasdataproperty(io::CSeis, name)

Return true if the data property called `name` exists in `io`.
"""
function hasdataproperty(io::CSeis, nm)
    for datapropkey in keys(io.dataproperties)
        if datapropkey == Symbol(nm)
            return true
        end
    end
    false
end

"""
    copy!(ioout::CSeis, hdrsout::AbstractMatrix{UInt8}, ioin::CSeis, hdrsin::AbstractMatrix{UInt8})

copy headers from `hdrsin` to `hdrsout`, and where `hdrsin` are headers from a frame of `ioin`, and
where `hdrsout` are headers from a frame of `ioout`.

# Example
```julia
ioin = csopen(AzContainer("mydataset-in-cs";storageaccount="mystorageaccount"))
ioout = csopen(AzContainer("mydataset-out-cs";storageaccount="mystorageaccount"))

hdrsin = readframehdrs(ioin, 1)
hdrsout = allocframehdrs(ioout)
copy!(ioout, hdrsout, ioin, hdrsin)
```
"""
function Base.copy!(ioout::CSeis, hdrsout::AbstractArray{UInt8,2}, ioin::CSeis, hdrsin::AbstractArray{UInt8,2})
    for trcpropout in ioout.traceproperties
        for trcpropin in ioin.traceproperties
            if trcpropout.def.label == trcpropin.def.label
                for i = 1:size(hdrsout,2)
                    set!(trcpropout, hdrsout, i, get(trcpropin, hdrsin, i))
                end
                break
            end
        end
    end
end

export
LogicalIndices,
DataProperty,
Geometry,
TracePropertyDef,
allocframe,
allocframetrcs,
allocframehdrs,
cache,
cache!,
dataproperty,
description!,
domains,
fold,
fold!,
geometry,
getframe,
getframehdrs,
getframetrcs,
cscreate,
csopen,
hasdataproperty,
headerlength,
lincs,
lstarts,
pincs,
prop,
propdefs,
props,
pstarts,
readhdrs,
readhdrs!,
readframe,
readframe!,
readframehdrs,
readframehdrs!,
readframetrcs,
readframetrcs!,
readtrcs,
readtrcs!,
stockdatatype,
stockdomain,
stockprop,
stockunit,
set!,
tracetype,
units,
writeframe

if !isdefined(Base, :get_extension)
    include("../ext/CvxCompressExt.jl")
    include("../ext/ZfpExt.jl")
end

end
