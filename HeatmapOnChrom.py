import svgwrite

def read_data(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            chrom, start, end, value = line.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((int(start), int(end), float(value)))
    return data

def get_chromosome_length(data):
    return max(end for _, end, _ in data)

def create_chromosome_shape(x, y, width, height):
    # Create a rounded rectangle shape for chromosomes
    radius = height / 2
    path = svgwrite.path.Path(d=f"M {x+radius},{y}")
    path.push(f"L {x+width-radius},{y}")
    path.push(f"A {radius},{radius} 0 0 1 {x+width-radius},{y+height}")
    path.push(f"L {x+radius},{y+height}")
    path.push(f"A {radius},{radius} 0 0 1 {x+radius},{y}")
    path.push("Z")
    return path

def extract_numeric_part(chrom):
    # Extract numeric part from chromosome name; return a large number for non-numeric parts
    numeric_part = ''.join(filter(str.isdigit, chrom))
    return int(numeric_part) if numeric_part else float('inf')  # Treat non-numeric as high value

def create_visualization(genes_file, output_file):
    genes_data = read_data(genes_file)

    # Sort chromosomes by their numerical order or alphabetically if they contain no numbers
    chromosomes = sorted(genes_data.keys(), key=extract_numeric_part)

    # Set up SVG canvas dimensions
    svg_width = 1200
    chr_height = 60
    chr_spacing = 40
    margin = 50
    svg_height = len(chromosomes) * (chr_height + chr_spacing) + 2 * margin

    dwg = svgwrite.Drawing(output_file, size=(svg_width, svg_height))

    # Add white background
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'), fill='white'))

    # Calculate max chromosome length for normalization
    max_chr_len = max(get_chromosome_length(genes_data[chrom]) for chrom in chromosomes)
    max_gene_value = max(value for chrom in genes_data for _, _, value in genes_data[chrom])

    # Define a reversed color scale (dark for max, light for min)
    colors = ['#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061']

    # Draw each chromosome with corresponding heatmap
    for i, chrom in enumerate(chromosomes):
        y = margin + i * (chr_height + chr_spacing)
        # Chromosome width adjusted to match the heatmap length exactly
        chr_width = ((svg_width - 2 * margin) * get_chromosome_length(genes_data[chrom])) / max_chr_len
        
        # Create and add the chromosome shape
        chr_shape = create_chromosome_shape(margin+33, y, chr_width, chr_height)
        chr_shape.fill('none').stroke('black', width=3)
        dwg.add(chr_shape)

        # Create a clip path to contain the heatmap within the chromosome shape
        clip_path = dwg.defs.add(dwg.clipPath(id=f'clip-{chrom}'))
        clip_path.add(chr_shape)

        # Add heatmap data within the chromosome shape
        heatmap_group = dwg.g(clip_path=f'url(#clip-{chrom})')
        for start, end, value in genes_data[chrom]:
            x = margin +33 + (start / get_chromosome_length(genes_data[chrom])) * chr_width
            width = max(1, ((end - start) / get_chromosome_length(genes_data[chrom])) * chr_width)
            color_index = min(int(value / max_gene_value * len(colors)), len(colors) - 1)
            heatmap_group.add(dwg.rect(insert=(x, y), size=(width, chr_height),
                                       fill=colors[color_index], stroke='none'))
        dwg.add(heatmap_group)

        # Add chromosome label close to the chromosome shape with increased font size
        dwg.add(dwg.text(chrom, insert=(margin + 25, y + chr_height/2 + 5),  # Position adjusted for closeness
                         font_size=30, fill='black', text_anchor='end'))  # Increased font size to 20

    # Add scale bar with further increased font size
    scale_y = svg_height - margin/2
    for i in range(6):
        x = margin +33 + i * (svg_width - 2 * margin) / 5
        dwg.add(dwg.line(start=(x, scale_y), end=(x, scale_y + 5), stroke='black', stroke_width=1))
        dwg.add(dwg.text(f'{i * max_chr_len / 5 / 1e6:.1f} Mb', insert=(x, scale_y + 20),
                         font_size=26, text_anchor='middle'))  # Increased font size to 18

    # Add heatmap scale with reversed colors and further increased font size
    heatmap_scale_width = 100
    heatmap_scale_height = 20
    for i, color in enumerate(colors):  # Colors from dark to light
        x = svg_width - margin - heatmap_scale_width + (i * heatmap_scale_width / len(colors))
        dwg.add(dwg.rect(insert=(x, margin/2), size=(heatmap_scale_width/len(colors), heatmap_scale_height),
                         fill=color, stroke='none'))
    dwg.add(dwg.text('Min', insert=(svg_width - margin - heatmap_scale_width, margin/2 + heatmap_scale_height + 15),
                     font_size=26, text_anchor='start'))  # Increased font size to 18
    dwg.add(dwg.text('Max', insert=(svg_width - margin, margin/2 + heatmap_scale_height + 15),
                     font_size=26, text_anchor='end'))  # Increased font size to 18

    dwg.save()

# Usage
create_visualization('AllChr2.cir', 'genome_visualization2.svg')

