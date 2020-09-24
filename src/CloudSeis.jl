module CloudSeis

using AbstractStorage, Blosc, Distributed, JSON, Random, TeaSeis

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

Base.Dict(dataproperty::DataProperty) = Dict("label"=>dataproperty.label,"value"=>dataproperty.value)
DataProperty(dataproperty::Dict) = DataProperty(dataproperty["label"], dataproperty["value"])

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
    Dict("name"=>extent.name, "container"=>JSON.parse(json(scrubsession(extent.container))), "firstframe"=>extent.frameindices[1], "lastframe"=>extent.frameindices[end])
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

# compression algorithms <--
abstract type AbstractCompressor end

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
        "library_version" => "0.7.0",
        "library_options" => Dict(
            "algorithm" => c.algorithm,
            "level" => c.level,
            "shuffle" => c.shuffle))
end
Base.copy(c::BloscCompressor) = BloscCompressor(c.algorithm, c.level, c.shuffle)

struct NotACompressor <: AbstractCompressor end
NotACompressor(d::Dict) = NotACompressor()
Dict(c::NotACompressor) = Dict("method" => "none")
Base.copy(c::NotACompressor) = NotACompressor()

function Compressor(d::Dict)
    method = get(d, "method", "") # enables backwards compatability
    if method == "blosc"
        return BloscCompressor(d)
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
* `compressor="none"` Compress the cache before writing to disk.  This is particularly useful for data with variable fold.  chooose from: ("none", "blosc")

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
        compressor = "none")
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
        compressor = compressor)
    if mode == "r" || mode == "r+"
        _kwargs = process_kwargs(;kwargs...)
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

`containers` is either of type `Container` or `Vector{<:Container}`.  Int
the former case, all extents are stored in a single container, and in the later case, extents
are shraded accross multiple containers.
"""
cscreate(containers::Union{Container,Vector{<:Container}}; kwargs...) = close(csopen(containers, "w"; kwargs...))

function csopen_read(containers::Vector{<:Container}, mode)
    description = JSON.parse(read(containers[1], "description.json", String))

    traceproperties = get_trace_properties(description)
    csopen_from_description(containers, mode, description, traceproperties)
end

function make_extents(containers::Vector{C}, nextents, foldername, nominal_frames_per_extent, remaining_frames) where {C<:Container}
    extents = Vector{Extent{C}}(undef, nextents)
    k = 1
    l = ceil(Int, log10(nextents))
    lastframe = 0
    for iextent = 1:nextents
        firstframe = lastframe + 1
        lastframe = firstframe + nominal_frames_per_extent
        if iextent > remaining_frames
            lastframe -= 1
        end
        extents[iextent] = Extent("$foldername/extent-$(lpad(iextent,l,'0'))", containers[k], firstframe:lastframe)
        k = k == length(containers) ? 1 : k + 1
    end
    extents
end

function csopen_write(containers::Vector{<:Container}, mode; kwargs...)
    ndim = length(kwargs[:axis_lengths])
    ndim == 0 && error("must specify axis_lengths")
    axis_propdefs = kwargs[:axis_propdefs]
    if length(axis_propdefs) == 0
        axis_propdefs = [stockprop[:SAMPLE], stockprop[:TRACE], stockprop[:FRAME], stockprop[:VOLUME], stockprop[:HYPRCUBE]][1:ndim]
        for i = 6:ndim
            push!(axis_propdefs, TracePropertyDef("DIM$i", "DIM$i", Int32, 1))
        end
    end
    axis_units = length(kwargs[:axis_units]) == 0 ? ["unknown" for i=1:ndim] : kwargs[:axis_units]
    axis_domains = length(kwargs[:axis_domains]) == 0 ? ["unkown" for i=1:ndim] : kwargs[:axis_domains]
    axis_pstarts = length(kwargs[:axis_pstarts]) == 0 ? [0.0 for i=1:ndim] : kwargs[:axis_pstarts]
    axis_pincs = length(kwargs[:axis_pincs]) == 0 ? [1.0 for i=1:ndim] : kwargs[:axis_pincs]

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
            bytes_total = prod(axis_lengths[2:end])*(axis_lengths[1]*sizeof(traceformat) + headerlength(traceproperties)) + sizeof(Int)*nframes
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
            "axis_pincs" => axis_pincs),
        "traceproperties"=>[Dict(traceproperty) for traceproperty in traceproperties],
        "dataproperties"=>[Dict(dataproperty) for dataproperty in kwargs[:dataproperties]],
        "extents"=>Dict.(extents),
        "compressor"=>Dict(kwargs[:compressor]))

    if kwargs[:geometry] != nothing
        merge(description, Dict("geometry"=>Dict(kwargs[:geometry])))
    end
    if length(kwargs[:dataproperties]) > 0
        merge(description, Dict("dataproperties"=>[Dict(dataproperty) for dataproperty in kwargs[:dataproperties]]))
    end

    mkpath.(containers)
    write(containers[1], "description.json", json(description, 1))
    io = csopen_from_description(containers, mode, description, traceproperties)
    empty!(io) # agressive emptying (over-writing) of existing data-set with same name.
    io
end

function process_kwargs(;kwargs...)
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
        compressor = (kwargs[:compressor] == "none" || kwargs[:compressor] == "") ? NotACompressor() : BloscCompressor("blosclz", 5, true)
    )
end

function process_kwargs_similarto(;kwargs...)
    io = csopen(kwargs[:similarto])

    local extents
    if kwargs[:frames_per_extent] == 0 && kwargs[:mbytes_per_extent] < 0 && length(kwargs[:extents]) == 0
        extents = [io.extents[i].frameindices for i=1:length(io.extents)]
    else
        extents = kwargs[:extents]
    end

    names = fieldnames(typeof(io.dataproperties))
    dataproperties = [DataProperty(String(names[i]), io.dataproperties[i]) for i=1:length(io.dataproperties)]

    (
        similarto = kwargs[:similarto],
        datatype = kwargs[:datatype] == "" ? io.datatype : kwargs[:datatype],
        force = kwargs[:force],
        traceformat = kwargs[:traceformat] == nothing ? io.traceformat : kwargs[:traceformat],
        byteorder = kwargs[:byteorder] == "" ? io.byteorder : kwargs[:byteorder],
        extents = extents,
        frames_per_extent = kwargs[:frames_per_extent],
        mbytes_per_extent = kwargs[:mbytes_per_extent] < 0 ? 1024 : kwargs[:mbytes_per_extent],
        geometry = kwargs[:geometry] == nothing ? io.geometry : kwargs[:geometry],
        tracepropertydefs = isempty(kwargs[:tracepropertydefs]) ? [io.traceproperties[i].def for i=1:length(io.traceproperties)] : kwargs[:tracepropertydefs],
        dataproperties = isempty(kwargs[:dataproperties]) ? dataproperties : kwargs[:dataproperties],
        axis_propdefs = isempty(kwargs[:axis_propdefs]) ? [io.axis_propdefs[i] for i=1:length(io.axis_propdefs)] : kwargs[:axis_propdefs],
        axis_units = isempty(kwargs[:axis_units]) ? [io.axis_units[i] for i=1:length(io.axis_units)] : kwargs[:axis_units],
        axis_domains = isempty(kwargs[:axis_domains]) ? [io.axis_domains[i] for i=1:length(io.axis_domains)] : kwargs[:axis_domains],
        axis_lengths = isempty(kwargs[:axis_lengths]) ? [io.axis_lengths[i] for i=1:length(io.axis_lengths)] : kwargs[:axis_lengths],
        axis_pstarts = isempty(kwargs[:axis_pstarts]) ? [io.axis_pstarts[i] for i=1:length(io.axis_pstarts)] : kwargs[:axis_pstarts],
        axis_pincs = isempty(kwargs[:axis_pincs]) ? [io.axis_pincs[i] for i=1:length(io.axis_pincs)] : kwargs[:axis_pincs],
        compressor = isempty(kwargs[:compressor]) ? similar(io.cache.compressor) : (kwargs[:compressor] == "none" ? NotACompressor() : BloscCompressor("blosclz", 5, true))
    )
end

function csopen_from_description(containers, mode, description, traceproperties)
    ndim = length(description["fileproperties"]["axis_lengths"])
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
        get_axis_propdefs(description, traceproperties),
        traceproperties,
        get_data_properties(description),
        get_geometry(description),
        [Extent(extent, containers) for extent in description["extents"]],
        Cache(Compressor(get(description, "compressor", Dict()))),
        headerlength(traceproperties))
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

function headerlength(traceproperties::NamedTuple)
    hdrlength = 0
    for traceproperty in traceproperties
        hdrlength += sizeof(eltype(traceproperty.def.format))*traceproperty.def.elementcount
    end
    hdrlength
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
    names = ntuple(i->Symbol(c[i]["label"]), length(c))
    values = [c[i]["value"] for i=1:length(c)]
    NamedTuple{names}(values)
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

function cache_from_file!(io::CSeis{T,N,NotACompressor}, extentindex) where {T,N}
    @debug "reading extent $extentindex from block-storage..."
    t = @elapsed begin
        io.cache.data = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end
    mb = length(io.cache.data)/1_000_000
    mbps = mb/t
    @debug "...data read ($mbps MB/s -- $mb MB)"
    nothing
end

function cache_from_file!(io::CSeis{T,N,BloscCompressor}, extentindex) where {T,N}
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

function cache!(io::CSeis, extentindex::Integer, force=false)
    if extentindex == io.cache.extentindex && io.cache.type == CACHE_ALL && !force
        return extentindex
    end

    if io.mode == "w" || io.mode == "r+"
        flush(io)
    end

    if isfile(io.extents[extentindex].container, io.extents[extentindex].name)
        cache_from_file!(io, extentindex)
    else
        @debug "creating cache..."
        io.cache.data = zeros(UInt8, cachesize(io, extentindex))
        @debug "...done creating cache."
    end
    io.cache.extentindex = extentindex
    io.cache.type = CACHE_ALL

    extentindex
end

function cache!(io::CSeis, idx::CartesianIndex, force=false)
    extentindex = extentindex_from_frameindex(io, idx)
    cache!(io, extentindex, force)
end

cache(io::CSeis) = io.cache

function cache_foldmap!(io::CSeis{T,N,NotACompressor}, extentindex::Integer, force=false) where {T,N}
    _partialcache_error() = error("partial caching is only allowed in read-only \"r\" mode")
    io.mode == "r" || _partialcache_error()

    if extentindex == io.cache.extentindex && io.cache.type ∈ (CACHE_FOLDMAP,CACHE_ALL) && !force
        return extentindex
    end

    @debug "reading foldmap for extent $extentindex from block-storage..."
    t = @elapsed begin
        io.cache.data = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, cachesize_foldmap(io, extentindex)))
    end
    mb = length(io.cache.data)/1_000_000
    mbps = mb/t
    @debug "...foldmap read ($mbps MB/s -- $mb MB)"
    io.cache.extentindex = extentindex
    io.cache.type = CACHE_FOLDMAP

    extentindex
end

function cache_foldmap!(io::CSeis{T,N,BloscCompressor}, extentindex::Integer, force=false) where {T,N}
    if extentindex == io.cache.extentindex && io.cache.type ∈ CACHE_ALL && !force
        return extentindex
    end
    # TODO - we can't to a partial cache here because of how the compression works.
    cache!(io, extentindex, force)
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

    function flush_blosc_buffer_byterange(iblock, block_size, block_remainder)
        isremainder = iblock <= block_remainder
        firstbyte = (iblock - 1)*block_size + (isremainder ? iblock : block_remainder + 1)
        lastbyte = firstbyte + (isremainder ? block_size : block_size - 1)
        firstbyte,lastbyte
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        Blosc.set_num_threads(Sys.CPU_THREADS)
        Blosc.set_compressor(io.cache.compressor.algorithm)

        maxbuffersize = 2_000_000_000
        cachesize = length(io.cache.data)
        nbuffers,nremainder = divrem(cachesize, maxbuffersize)
        nremainder > 0 && (nbuffers += 1)
        buffersize,nremainder = divrem(cachesize, nbuffers)

        buffer_lengths = zeros(Int, nbuffers)
        cdata = zeros(UInt8, 8*(1 + nbuffers)) # store number of buffers and length of each compressed buffer
        cdata_io = IOBuffer(cdata; read=true, write=true)
        write(cdata_io, nbuffers)
        close(cdata_io)
        for ibuffer = 1:nbuffers
            firstbyte,lastbyte = flush_blosc_buffer_byterange(ibuffer, buffersize, nremainder)
            _data = unsafe_wrap(Array, pointer(io.cache.data)+(firstbyte-1), (lastbyte-firstbyte+1,))
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

function foldmap(io::CSeis)
    nbytes = sizeof(Int)*length(io.extents[io.cache.extentindex].frameindices)
    reinterpret(Int, view(io.cache.data, 1:nbytes))
end

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
    fmap[LinearIndices(io)[idx] - io.extents[extentindex].frameindices[1] + 1]
end
TeaSeis.fold(io::CSeis, idx...) = fold(io, CartesianIndex(idx))

function TeaSeis.fold(io::CSeis, hdrs::AbstractArray{UInt8,2})
    fld = 0
    for i = 1:size(hdrs,2)
        if get(prop(io,stockprop[:TRC_TYPE]), hdrs, i) == tracetype[:live]
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
    fmap[LinearIndices(io)[idx] - io.extents[extentindex].frameindices[1] + 1] = fld
    nothing
end
fold!(io::CSeis, fld, idx...) = fold!(io, fld, CartesianIndex(idx))

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

function extentindex_from_frameindex(io, idx::CartesianIndex)
    frameidx = LinearIndices(io)[idx]

    # TODO - optimize me
    for (i,extent) in enumerate(io.extents)
        if frameidx ∈ extent.frameindices
            return i
        end
    end
    error("can't find extent")
end

function trcsoffset(io::CSeis, extentindex, idx)
    nsamples = size(io,1)
    ntraces = size(io,2)
    nframes = length(io.extents[extentindex].frameindices)
    frstframe = io.extents[extentindex].frameindices[1]

    frameidx = LinearIndices(io)[idx]
    sizeof(Int)*nframes + io.hdrlength*ntraces*nframes + sizeof(io.traceformat)*nsamples*ntraces*(frameidx-frstframe) + 1
end

function hdrsoffset(io::CSeis, extentindex, idx)
    ntraces = size(io,2)
    nframes = length(io.extents[extentindex].frameindices)
    frstframe = io.extents[extentindex].frameindices[1]

    frameidx = LinearIndices(io)[idx]
    sizeof(Int)*nframes + io.hdrlength*ntraces*(frameidx-frstframe) + 1
end

function getframetrcs(io::CSeis, idx::CartesianIndex)
    extentindex = cache!(io, idx)
    data = io.cache.data

    frstbyte = trcsoffset(io, extentindex, idx)
    lastbyte = frstbyte + size(io,2)*size(io,1)*sizeof(io.traceformat) - 1
    reshape(reinterpret(io.traceformat, view(data,frstbyte:lastbyte)), :, size(io,2))
end
getframetrcs(io, idx...) = getframetrcs(io, CartesianIndex(idx))

function getframehdrs(io::CSeis, idx::CartesianIndex)
    extentindex = cache!(io, idx)
    data = io.cache.data

    frstbyte = hdrsoffset(io, extentindex, idx)
    lastbyte = frstbyte + size(io,2)*io.hdrlength - 1
    reshape(view(data,frstbyte:lastbyte), :, size(io,2))
end
getframehdrs(io::CSeis, idx...) = getframehdrs(io, CartesianIndex(idx))

getframe(io::CSeis, idx) = getframetrcs(io, idx), getframehdrs(io, idx)

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

"""
    readframetrcs!(io::CSeis, trcs, idx...)

Read traces from `io` into `trcs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.
"""
function TeaSeis.readframetrcs!(io::CSeis, trcs::AbstractArray, idx::CartesianIndex)
    _trcs = getframetrcs(io, idx)
    fld = fold(io, idx)
    copyto!(trcs, 1, _trcs, 1, size(io,1)*fld)
end

"""
    trcs = readframetrcs(io::CSeis, idx...)

Read traces from `io` for the frame `idx...`.  `idx..` can either be
integer(s) or a `CartesianIndex`.
"""
TeaSeis.readframetrcs(io::CSeis, idx::CartesianIndex) = readframetrcs!(io, allocframetrcs(io), idx)

TeaSeis.readframetrcs!(io::CSeis, trcs::AbstractArray, idx...) = readframetrcs!(io, trcs, CartesianIndex(idx))
TeaSeis.readframetrcs(io::CSeis, idx...) = readframetrcs(io, CartesianIndex(idx))

"""
    readframehdrs!(io::CSeis, trcs, idx...)

Read headers from `io` into `hdrs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.
"""
function TeaSeis.readframehdrs!(io::CSeis, hdrs::AbstractArray{UInt8,2}, idx::CartesianIndex)
    _hdrs = getframehdrs(io, idx)
    fld = fold(io, idx)
    copyto!(hdrs, 1, _hdrs, 1, io.hdrlength*fld)
end

"""
    hdrs = readframehdrs(io::CSeis, idx...)

Read headers from `io` for the frame `idx...`.  `idx..` can either be
integer(s) or a `CartesianIndex`.
"""
TeaSeis.readframehdrs(io::CSeis, idx::CartesianIndex) = readframehdrs!(io, allocframehdrs(io), idx)

TeaSeis.readframehdrs!(io::CSeis, hdrs::AbstractArray{UInt8,2}, idx...) = readframehdrs!(io, hdrs, CartesianIndex(idx))
TeaSeis.readframehdrs(io::CSeis, idx...) = readframehdrs(io, CartesianIndex(idx))

"""
    readframe!(io::CSeis, trcs, hdrs, idx...)

Read traces and headers from `io` into `trcs::Matrix`, and `hdrs::Matrix` for the frame `idx...`.
`idx...` can either be integer(s) or a `CartesianIndex`.
"""
TeaSeis.readframe!(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2}, idx::CartesianIndex) = readframetrcs!(io, trcs, idx), readframehdrs!(io, hdrs, idx)

"""
    trcs,hdrs = readframe(io::CSeis, trcs, hdrs, idx...)

Read traces and headers from `io` for the frame `idx...`.  `idx...` can
either be integer(s) or a `CartesianIndex`.
"""
TeaSeis.readframe(io::CSeis, idx::CartesianIndex) = readframe!(io, allocframetrcs(io), allocframehdrs(io), idx)
TeaSeis.readframe!(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2}, idx...) = readframe!(io, trcs, hdrs, CartesianIndex(idx))
TeaSeis.readframe(io::CSeis, idx...) = readframe(io, CartesianIndex(idx))

"""
    writeframe(io::CSeis, trcs, hdrs)

Write a frame to `io`.  The location of the frame is determined
by the axis headers set in `hdrs`.
"""
function TeaSeis.writeframe(io::CSeis, trcs::AbstractArray, hdrs::AbstractArray{UInt8,2})
    idx = CartesianIndex(ntuple(i->get(prop(io, io.axis_propdefs[i+2]), hdrs, 1), ndims(io)-2))
    cache!(io, idx)
    fld = fold(io, hdrs)
    fold!(io, fld, idx)
    _trcs = getframetrcs(io, idx)
    _hdrs = getframehdrs(io,idx)
    copyto!(_trcs, 1, trcs, 1, fld*size(io,1))
    copyto!(_hdrs, 1, hdrs, 1, fld*io.hdrlength)
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
    _trcs,hdrs = getframe(io, idx)
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
        idx = parseindex(rng, idx_n)
        frmtrcs = getframetrcs(io, idx)
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
trcs = readtrcs!(io, 1:50, :, 1)
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
        idx = parseindex(rng, idx_n)
        frmhdrs = getframehdrs(io, idx)
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
        idx = parseindex(_rng, idx_n)

        frmtrcs = getframetrcs(io, idx)
        for (itrc,trc) in enumerate(trcrng), (ismp,smp) in enumerate(smprng)
            frmtrcs[smp,trc] = trcs[ismp,itrc,idx_n]
        end

        for itrace = 1:io.axis_lengths[2]
            set!(props[2], frmhdrs, itrace, itrace)
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
* frames_per_extent=0` set the size of each output extent via the number of frames allowed in each extent.
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

"""
    cp(src::CSeis, dst::Container)

copy a CloudSeis data-set to `dst`.
"""
function Base.cp(src::CSeis, dst::Container)
    iodst = csopen(dst, "w", similarto=dirname(src.containers[1]))
    trcs, hdrs = allocframe(iodst)
    for idx in CartesianIndices(size(io)[3:end])
        fld = readframe!(src, trcs, hdrs, idx)
        if fld > 0
            writeframe(iodst, trcs, hdrs)
        end
    end
    close(iodst)
end

"""
    mv(src::CSeis, dst::Container)

move a CloudSeis data-set to `dst`.
"""
function Base.mv(src::CSeis, dst::Container)
    cp(src, dst, project=project, location=location)
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
    dataproperty(io::CSeis, name)

Return the CloudSeis data property with `name`
that is in the CloudSeis data-set `io`.
"""
function dataproperty(io::CSeis, nm)
    try
        return io.dataproperties[Symbol(nm)]
    catch
        error("data property with label $nm does not exist")
    end
end

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
DataProperty,
Geometry,
TracePropertyDef,
allocframe,
allocframetrcs,
allocframehdrs,
cache,
cache!,
dataproperty,
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

end
