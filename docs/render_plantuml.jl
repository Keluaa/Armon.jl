
using CodecZlib
using HTTP


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


function encode_to_plant_uml_base64(filename::String)
    open(filename, "r") do file
        return encode_to_plant_uml_base64(file) 
    end
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
    return
end


function render_plantuml_files(files::Vector{String})
    @info "Rendering PlantUML files to SVG"
    svg_files = map(files) do file
        file_path = joinpath(@__DIR__, file)
        encoded_file = encode_to_plant_uml_base64(file_path)
        svg_file = first(splitext(file_path)) * ".svg"
        send_to_render_server(encoded_file, svg_file)
    end
    return svg_files
end
