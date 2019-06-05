from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
import os
directory = '/gpfs/ts0/home/bt273/BIOS-SCOPE/metag/ashley/genemaps'

for filename in os.listdir(directory):
    if filename.endswith(".gbk"):
        new_name=os.path.splitext(filename)[0]
        print(filename)
        record = next(SeqIO.parse(filename, "genbank"))
        gd_diagram = GenomeDiagram.Diagram(new_name)
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()
        
        for feature in record.features:
            if feature.type != "gene":
                continue
            if len(gd_feature_set) % 2 == 0:
                color = colors.blue
            else:
                color = colors.lightblue
            gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_size = 10)
            gd_diagram.draw(format="linear", orientation="landscape", pagesize='A1',fragments=25, start=0, end=len(record))
            gd_diagram.write(new_name+".pdf", "PDF")
