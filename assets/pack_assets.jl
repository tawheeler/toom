using Images
using FileIO

const N_BYTES_TOC_NAME = 16

wall_textures = load("textures.tif")
floor_textures = load("floor_textures.tif")

# The output format is:
# HEADER:
#    "TOOM"   - 4 chars
# Data:
#   necessary data, laid out as needed.
# Table of Contents:
#   array of table of content entries:
#      u32       offset in file # number of bytes past 'TOOM' to read at. First entry will have offset 0
#      char[16]  name           # null-terminated string label, e.g. "floor_textures"
#   u32 n_toc_entries = number of table of content entries

struct TableOfContentsEntry
	name::String
	offset_in_file::UInt32
end

function write_image(img, output::IOStream; column_major::Bool=true)::UInt32
	n_bytes_written = zero(UInt32)

	# Write the number of pixels as a UInt32
	n_pixels = UInt32(length(img))
	write(output, n_pixels::UInt32)
	n_bytes_written += 4

	# Write the number of pixels per column as a UInt32
	n_pix_per_column = UInt32(size(img)[1])
	write(output, n_pix_per_column::UInt32)
	n_bytes_written += 4

	# Write whether the pixels are in column major order
	write(output, UInt8(column_major))
	n_bytes_written += 1

	# Write the image top to bottom and left to right, as ABGR little endian
	if column_major
		# a11, a21, a31, a12, a22, a32, a31, a23, a33
		for col in 1:size(img)[2]
			for row in 1:n_pix_per_column
				write(output, (img[row, col].r.i)::UInt8)
				write(output, (img[row, col].g.i)::UInt8)
				write(output, (img[row, col].b.i)::UInt8)
				write(output, 0xFF);
				n_bytes_written += 4
			end
		end
	else
		# a11, a12, a13, a21, a22, a23, a31, a32, a33
		for row in 1:n_pix_per_column
			for col in 1:size(img)[2]
				write(output, (img[row, col].r.i)::UInt8)
				write(output, (img[row, col].g.i)::UInt8)
				write(output, (img[row, col].b.i)::UInt8)
				write(output, 0xFF);
				n_bytes_written += 4
			end
		end
	end

	return n_bytes_written
end

offset_in_file = zero(UInt32)
table_of_contents_entries = TableOfContentsEntry[]

output_file = "assets.bin"
open(output_file, "w") do output
	write(output, 'T')
	write(output, 'O')
	write(output, 'O')
	write(output, 'M')

	global offset_in_file
	global table_of_contents_entries

	# ADD CONTENT ---------------------------------------
	push!(table_of_contents_entries, TableOfContentsEntry(
		"wall_textures",
		offset_in_file
	))
	offset_in_file += write_image(wall_textures, output, column_major=true)

	push!(table_of_contents_entries, TableOfContentsEntry(
		"floor_textures",
		offset_in_file
	))
	offset_in_file += write_image(floor_textures, output, column_major=true)

	# WRITE TABLE OF CONTENTS -----------------------------

	for entry in table_of_contents_entries
		write(output, entry.offset_in_file)
		offset_in_file += 4

		for i in 1:N_BYTES_TOC_NAME-1
			c = '\0'
			if i < length(entry.name)
				c = entry.name[i]
			end
			write(output, c)
		end
		write(output, '\0') # Always null-terminate
		offset_in_file += N_BYTES_TOC_NAME
	end

	write(output, UInt32(length(table_of_contents_entries)))
	offset_in_file += 4
end

println("Done!")
println("Total bytecount: ", offset_in_file + 4)