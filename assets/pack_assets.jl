using Images
using FileIO

img = load("textures.tif")

output_file = "assets.bin"
open(output_file, "w") do output
	write(output, 'T')
	write(output, 'O')
	write(output, 'O')
	write(output, 'M')

	# Write the number of pixels as a UInt32
	n_pixels = UInt32(length(img))
	write(output, n_pixels::UInt32)

	# Write the number of pixels per column as a UInt32
	n_pix_per_column = UInt32(size(img)[1])
	write(output, n_pix_per_column::UInt32)

	# Write the image top to bottom and left to right, in BGR order (no transparency)
	for col in 1:size(img)[2]
		for row in 1:n_pix_per_column
			write(output, (img[row, col].r.i)::UInt8)
			write(output, (img[row, col].g.i)::UInt8)
			write(output, (img[row, col].b.i)::UInt8)
			write(output, 0xFF);
		end
	end
end

println("Done!")