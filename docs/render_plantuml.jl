module DocumenterPlantUML

using CodecZlib
using HTTP
using Documenter

export PlantUML


# plantUML uses a custom version of base64, for historical reasons
# See https://plantuml.com/en/text-encoding
const plantUML_base64 = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-_"
encode_b64(b::UInt8) = @inbounds plantUML_base64[(b & 0x3F) + 1]


function encode_b64(bytes::Vector{UInt8})
    encoded_bytes = Vector{UInt8}(undef, ceil(Int, length(bytes) * 4/3))
    next_pos = 1

    len = length(bytes)
    for i in 1:3:lastindex(bytes)
        b1 = bytes[i]
        b2 = i+1 ≤ len ? bytes[i+1] : 0x00
        b3 = i+2 ≤ len ? bytes[i+2] : 0x00

        # 3 bytes to 4 chars
        c1 = encode_b64((b1 & 0b1111_1100) >> 2)
        c2 = encode_b64(((b1 & 0b0000_0011) << 4) | ((b2 & 0b1111_0000) >> 4))
        c3 = encode_b64(((b2 & 0b0000_1111) << 2) | ((b3 & 0b1100_0000) >> 6))
        c4 = encode_b64(b3 & 0b0011_1111)

        encoded_bytes[next_pos] = c1; next_pos += 1
        encoded_bytes[next_pos] = c2; next_pos += 1
        i+1 ≤ len && (encoded_bytes[next_pos] = c3; next_pos += 1)
        i+2 ≤ len && (encoded_bytes[next_pos] = c4; next_pos += 1)
    end

    # no padding, I don't know if that is supported by the server
    return String(@view encoded_bytes[1:next_pos-1])
end


function encode_to_plant_uml_base64(io::IO)
    deflated = read(DeflateCompressorStream(io))
    in_base64 = encode_b64(deflated)
    return in_base64
end


function send_to_render_server(encoded_file::String, svg_filepath)
    if length(encoded_file) > 4000
        # Avoid unhelpful HTTP request errors when the file is too large
        error("PUML file is too large, encoded length is: $(length(encoded_file)), limit is ~4000")
    end

    # See https://plantuml.com/en/server
    HTTP.open(:GET, "https://plantuml.com/plantuml/svg/" * encoded_file) do http
        open(svg_filepath, "w") do svg_file
            write(svg_file, http)
        end
    end

    return svg_filepath
end


function render_puml(puml_source::IO, svg_path)
    encoded_file = encode_to_plant_uml_base64(puml_source)
    send_to_render_server(encoded_file, svg_path)
end

render_puml(puml_source::AbstractString, svg_path) = render_puml(IOBuffer(puml_source), svg_path)


struct PlantUML <: Documenter.Plugin
    # source '.puml' path to (modified PUML source, '.svg' result path)
    diagrams::Dict{String, Tuple{String, String}}

    PlantUML() = new(Dict())
end


abstract type PlantUMLBlocks <: Documenter.Expanders.NestedExpanderPipeline end

Documenter.Selectors.order(::Type{PlantUMLBlocks}) = 11.1  # After @raw blocks
Documenter.Selectors.matcher(::Type{PlantUMLBlocks}, node, page, doc) = Documenter.iscode(node, r"^@puml")

function Documenter.Selectors.runner(::Type{PlantUMLBlocks}, node, page, doc)
    @assert node.element isa Documenter.MarkdownAST.CodeBlock
    x = node.element

    m = match(r"@puml( .+)?$", x.info)
    m === nothing && error("invalid '@puml [PlantUML file path]' syntax: $(x.info)")

    plant_uml = Documenter.getplugin(doc, PlantUML)

    if !isnothing(m[1])
        file_path = joinpath(doc.user.source, strip(m[1]))
        !isfile(file_path) && @error "file '$file_path' does not exist"
        puml_source = read(file_path, String)
    else
        # idea: inline PlantUML script in a documenter file
        # simply put the `x.code` in a tmp file and render it
        @error "missing PlantUML source file path"
    end

    svg_file = splitext(splitdir(file_path)[end])[1] * ".svg"
    svg_path = joinpath(doc.user.build, "assets", svg_file)
    rel_svg_path = relpath(svg_path, doc.user.build)
    plant_uml.diagrams[file_path] = (svg_path, puml_source)

    # Placing our SVG result in a `embed` tag allow selectable text and clickable links
    # without extra effort
    code = """
    <embed src="$rel_svg_path" />
    """

    node.element = Documenter.RawNode(:html, code)
end


abstract type PlantUMLRender <: Documenter.Builder.DocumentPipeline end

Documenter.Selectors.order(::Type{PlantUMLRender}) = 6.1  # After rendering the document

function Documenter.Selectors.runner(::Type{PlantUMLRender}, doc::Documenter.Document)
    Documenter.is_doctest_only(doc, "PlantUMLRender") && return
    plant_uml = Documenter.getplugin(doc, PlantUML)
    @info "PlantUML: Rendering $(length(plant_uml.diagrams)) diagrams"
    for (puml_file, (target_svg, puml_source)) in plant_uml.diagrams
        try
            render_puml(puml_source, target_svg)
        catch e
            @warn "Failed to render PlantUML file: '$puml_file'"
            rethrow(e)
        end
    end
end

end
