module DocumenterPlantUML

using CodecZlib
using HTTP
using Documenter
using MarkdownAST

import Documenter: HTMLWriter
import MarkdownAST: Node

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


function send_to_render_server(encoded_file::String, svg_filepath, dark_mode=false)
    if length(encoded_file) > 4000
        # Avoid unhelpful HTTP request errors when the file is too large
        error("PUML file is too large, encoded length is: $(length(encoded_file)), limit is ~4000")
    end

    # See https://plantuml.com/en/server
    d = dark_mode ? "d" : ""
    HTTP.open(:GET, "https://plantuml.com/plantuml/$(d)svg/" * encoded_file) do http
        open(svg_filepath, "w") do svg_file
            write(svg_file, http)
        end
    end

    return svg_filepath
end


function render_puml(puml_source::IO, svg_path, dark_mode=false)
    encoded_file = encode_to_plant_uml_base64(puml_source)
    send_to_render_server(encoded_file, svg_path, dark_mode)
end

render_puml(puml_source::AbstractString, svg_path, dark_mode=false) =
    render_puml(IOBuffer(puml_source), svg_path, dark_mode)


"""
    PlantUML(; no_render=false, silent=false, no_links=false, short_links=true, dark_mode=true)

Documenter plugin to control how PlantUML diagrams are rendered.

Options:
 - `no_render`: skip rendering of diagrams
 - `silent`: print no messages when rendering
 - `no_links`: don't replace links in the PlantUML source
 - `short_links`: allow replacing links of the form `[[\`SomeJuliaObj\`]]`
 - `dark_mode`: render diagrams with a dark background
"""
Base.@kwdef struct PlantUML <: Documenter.Plugin
    no_render   :: Bool = false
    silent      :: Bool = false
    no_links    :: Bool = false
    short_links :: Bool = true
    dark_mode   :: Bool = true
end


mutable struct PlantUMLDiagram
    path        :: String  # source '.puml' path
    svg_file    :: String  # '.svg' result path
    src_md_file :: String  # '.md' file including the diagram
    puml_source :: String  # PUML source
    height      :: String  # HTML height attribute. Special value of "exact" will match the SVG height.
    width       :: String  # HTML width attribute. Special value of "exact" will match the SVG width.
end


struct PlantUMLDiagramBlock <: MarkdownAST.AbstractBlock
    diagram::PlantUMLDiagram
end

MarkdownAST.iscontainer(::PlantUMLDiagramBlock) = true


function Documenter.MDFlatten.mdflatten(io, node::Node, e::PlantUMLDiagramBlock)
    print(io, "(PlantUML diagram ", e.diagram.path, " :")
    Documenter.MDFlatten.mdflatten(io, node.children)
    print(io, ")")
end


function HTMLWriter.domify(dctx::HTMLWriter.DCtx, node::Node, element::PlantUMLDiagramBlock)
    raw_html_diagram_embed = first(node.children)
    link_list = last(node.children)  # This dummy list is only here to delegate the resolving of URLs

    puml = Documenter.getplugin(dctx.ctx.doc, PlantUML)
    diagram = element.diagram
    if !puml.no_links
        # Get the URLs from the links embedded in the diagram 
        resolved_urls = map(link_list.children) do link_item
            resolved_link_node = first(first(link_item.children).children)  # item node -> paragraph node -> child
            resolved_link = resolved_link_node.element
            return HTMLWriter.filehref(dctx, node, resolved_link)
        end
        final_puml_source = replace_all_links(diagram.puml_source, resolved_urls, puml)
    else
        final_puml_source = diagram.puml_source
    end

    svg_path = render_puml(dctx.ctx, diagram.path, diagram.svg_file, diagram.src_md_file, final_puml_source, puml.dark_mode)

    if diagram.height == diagram.width == "exact"
        scrollbar_size = "0px"  # both use the exact dimension, no scrollbar is needed
    else
        scrollbar_size = "17px"
    end

    height_attr = diagram.height == "exact" ? find_svg_attribute(svg_path, "height") : diagram.height
    width_attr  = diagram.width  == "exact" ? find_svg_attribute(svg_path, "width")  : diagram.width

    # Placing our SVG result in an `iframe` allow selectable text and clickable links
    # without extra effort. `embed` would also work but also prevent scrollbars.
    # The `iframe` is why we must parse the SVG source in order to find its dimensions
    # when using "exact" size, as it is not possible for attributes to use the size of the
    # embedded content. Parsing the SVG seems like the easiest way to solve this.
    iframe_html = """
    <iframe src="$(diagram.svg_file)"
        overflow-x="auto" overflow-y="auto"
        style="width: calc($width_attr + $scrollbar_size); height: calc($height_attr + $scrollbar_size);"
        type="text/svg">
    </iframe>
    """

    raw_html_diagram_embed.element = Documenter.RawNode(:html, iframe_html)

    # Only domify the `@raw html <iframe ... />` node
    return HTMLWriter.domify(dctx, raw_html_diagram_embed, raw_html_diagram_embed.element)
end


function find_svg_attribute(svg_path, attribute)
    attr_re = Regex("\\Q$attribute\\E=\"([^\"]+)\"")
    svg_source = read(svg_path, String)
    m = match(attr_re, svg_source)
    if isnothing(m)
        error("could not find '$attribute' attribute in svg source at $svg_path")
    end
    return m[1]
end


function render_puml(ctx, puml_file, target_svg, src_md_page, puml_source, dark_mode)
    page_dir = dirname(HTMLWriter.get_url(ctx, src_md_page))
    page_dir = joinpath(ctx.doc.user.build, page_dir)
    !isdir(page_dir) && mkpath(page_dir)  # We are here before `write_html` is called
    svg_path = joinpath(page_dir, target_svg)

    plant_uml = Documenter.getplugin(ctx.doc, PlantUML)
    !plant_uml.silent && @info "PlantUML: rendering '$puml_file'"
    !plant_uml.no_render && render_puml(puml_source, svg_path, dark_mode)

    return svg_path
end


abstract type PlantUMLBlocks <: Documenter.Expanders.NestedExpanderPipeline end

Documenter.Selectors.order(::Type{PlantUMLBlocks}) = 11.1  # After @raw blocks
Documenter.Selectors.matcher(::Type{PlantUMLBlocks}, node, page, doc) = Documenter.iscode(node, r"^@puml")

function Documenter.Selectors.runner(::Type{PlantUMLBlocks}, node, page, doc)
    @assert node.element isa MarkdownAST.CodeBlock
    x = node.element

    puml_plugin = Documenter.getplugin(doc, PlantUML)

    m = match(r"""@puml ("[^"]+"|\S+)((?> +[\w=]+)*)$""", x.info)
    m === nothing && error("invalid '@puml <PlantUML file path> [option=value]' syntax: $(x.info)")

    if !isnothing(m[1])
        file_path = joinpath(doc.user.source, strip(m[1]))
        !isfile(file_path) && @error "file '$file_path' does not exist"
        puml_source = read(file_path, String)
    else
        # idea: inline PlantUML script in a documenter file
        # simply put the `x.code` in a tmp file and render it
        @error "missing PlantUML source file path"
    end

    # The default dimensions of "100%"x"exact" will display the SVG in its whole height
    # with no scrollbar, while keeping the width of the rest of the documentation.
    height = "exact"
    width = "100%"
    options = split(m[2]; keepempty=false)
    for raw_option in options
        option, value = split(raw_option, "="; limit=2)
        if     option == "height"  height = value
        elseif option == "width"   width  = value
        else
            error("unknown `@puml` option: $raw_option")
        end
    end

    svg_file = splitext(splitdir(file_path)[end])[1] * ".svg"
    src_md_file = relpath(page.source, doc.user.source)

    diagram = PlantUMLDiagram(file_path, svg_file, src_md_file, puml_source, height, width)

    # Since the links are embedded in PlantUML, Documenter cannot find and resolve them alone.
    # To do this, we place the diagram in a temporary Markdown AST node, containing the
    # real raw HTML code AND a Markdown list of all links in the PlantUML source. The links
    # can then be resolved by Documenter automatically.
    node.element = PlantUMLDiagramBlock(diagram)
    push!(node.children, Node(Documenter.RawNode(:html, "")))  # code is set later
    !puml_plugin.no_links && push!(node.children, extract_all_links(puml_source, puml_plugin))
end


# To avoid any ambiguities, we only match links enclosed in quotes. This way the user
# can place links with any syntax (including Documenter's) while staying compatible with
# PlantUML's syntax and live renderers of diagrams
# Should match any link of the form `[["[julia](@ref link)"blabla]]`
# First capture group is `julia`, second is `@ref link`, third is the rest of the PlantUML
# link syntax, which shouldn't be touched.
# The fourth capture group is only used for other regex.
const PLANT_UML_DOCUMENTER_LINK_REGEX = r"""\[\["\[([^"]+)\]\(([^"]+)\)"([^\]]*)()\]\]"""
# Same as above but also matches `[[`julia`]]`, as shorthand for `[["[`julia`](@ref)" julia]]`
const PLANT_UML_SHORT_DOCUMENTER_LINK_REGEX = r"""\[\[(?>"\[([^"]+)\]\(([^"]+)\)"|`([^`]+)`)([^\]]*)\]\]"""

function link_regex(p::PlantUML)
    if p.short_links
        return PLANT_UML_SHORT_DOCUMENTER_LINK_REGEX
    else
        return PLANT_UML_DOCUMENTER_LINK_REGEX
    end
end


function extract_all_links(puml_source, puml::PlantUML)
    link_list = Node(MarkdownAST.List(:bullet, true))

    for link_match in eachmatch(link_regex(puml), puml_source)
        if puml.short_links && !isnothing(link_match[3])
            obj = link_match[3]
            ref = "@ref"
        else
            obj = link_match[1]
            ref = link_match[2]
        end

        # Build a ` - [link](@ref)` node
        item_node = MarkdownAST.@ast MarkdownAST.Item() do
            MarkdownAST.Paragraph() do
                MarkdownAST.Link(ref, "") do
                    MarkdownAST.Text(obj)
                end
            end
        end
        push!(link_list.children, item_node)
    end

    return link_list
end


function replace_all_links(puml_source, resolved_urls, puml::PlantUML)
    i = 1
    function substitute_link(matched_link)
        resolved_url = resolved_urls[i]
        i += 1
        m = match(link_regex(puml), matched_link)
        if puml.short_links
            if !isnothing(m[3]) && isempty(m[4])
                # [[`Obj`]] => [["link" Obj]]
                return """[["$resolved_url" $(m[3])]]"""  # the space is needed (see PlantUML link syntax)
            else
                # [[`Obj` blabla]] => [["link" blabla]]
                # OR
                # [["[Obj](@ref)" blabla]] => [["link" blabla]]
                return """[["$resolved_url"$(m[4])]]"""
            end
        else
            # [["[Obj](@ref)" blabla]] => [["link" blabla]]
            return """[["$resolved_url"$(m[3])]]"""
        end
    end

    return replace(puml_source, link_regex(puml) => substitute_link)
end

end
