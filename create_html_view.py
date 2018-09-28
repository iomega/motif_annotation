"""
Script to create html visualization of MAGMa annotated MS2LDA motifs
It reads a json file create by annotate_motis.py

Usage: python create_html_view.py <annotation.json >annotation.html
"""

import sys, json


def write_header():
    print '''
<html>
<head>
<title>MS2LDA_MAGMa</title>
<meta charset="utf-8"/>
<script type="text/javascript">    
'''
    print open('ChemDoodleWeb.js').read()
    print '''
</script>

<style type="text/css">

a.tooltip {
    position: relative; 
    top: 0px; 
    left: 0px;
}

a.tooltip:hover span {
    opacity: 1; 
    visibility: visible;
}

a.tooltip span {
    padding: 10px;
    top: 0px;
    left: 150px;     
    background-color: #EEEEEE; 
    color: #000000;
    height: auto;
    width: auto;
    border-radius: 5px; 
    opacity: 1; 
    position: absolute;
    visibility: hidden;
}
</style>

<script  type="text/javascript">
function show_mol(id, molblock, fragatoms="", type="loss") {
  var viewer = new ChemDoodle.ViewerCanvas(id, 200, 150);
  viewer.specs.bonds_width_2D = .6;
  viewer.specs.bonds_saturationWidthAbs_2D = 2.6;
  viewer.specs.bonds_hashSpacing_2D = 2.5;
  viewer.specs.atoms_font_size_2D = 10;
  viewer.specs.atoms_font_families_2D = ['Helvetica', 'Arial', 'sans-serif'];
  viewer.specs.atoms_displayTerminalCarbonLabels_2D = false;
  viewer.specs.atoms_implicitHydrogens_2D = false;
  var molecule = ChemDoodle.readMOL(molblock);
  molecule.scaleToAverageBondLength(14.4);
  var lossspecs = new ChemDoodle.structures.VisualSpecifications();
  lossspecs.bonds_width_2D = 2.0;
  if (type == "loss") {
      lossspecs.bonds_color = 'red';
      lossspecs.atoms_color = 'red';
      }
  else {
      lossspecs.bonds_color = 'cyan';
      lossspecs.atoms_color = 'cyan';
      }
  if (fragatoms != "") {
      var fragmentAtoms = fragatoms.split(',');
      molecule.bonds.forEach(function(b) {
        // only color bond black if both atoms are in fragmentAtoms
        a1 = molecule.atoms.indexOf(b.a1);
        a2 = molecule.atoms.indexOf(b.a2);
        if (a1 != -1 && fragmentAtoms.indexOf(a1 + '') == -1) {
          molecule.atoms[a1].specs = lossspecs;
          b.specs = lossspecs;
          }
        if (a2 != -1 && fragmentAtoms.indexOf(a2 + '') == -1) {
          molecule.atoms[a2].specs = lossspecs;
          b.specs = lossspecs;
          }
        });
      };
  viewer.loadMolecule(molecule)
  }
</script>
</head>
'''


def write_motif_annotations():
    viewid = 0
    for motif in motifannotations:
        print '<p style="clear: both"><br><b><u>{}</u><br>{}</b></p>'.format(motif['name'], '<br>'.join(motif['annotations']).encode('utf-8'))
        for feature in motif['features']:
            print '<p style="clear: both"><br><u>{}, {}</u><br></p>'.format(feature['name'], feature['probability'])
            for substr in feature['substr']:
                if feature['type'] == 'fragment' and 'mol' in substr:
                    viewid += 1
                    viewerid = 'view'+str(viewid) 
                    struct = '<script>show_mol("{}","{}")</script>'.format(viewerid, substr['mol'].replace('\n','\\n'))
                else:
                    struct = substr['smiles']
                print '<div style="clear: both"><button class="collapsible" style="float: left">{} ({})</button>'.format(struct, len(substr['matches']))
                print '<div class="content" style="display: none"><table style="background-color: #eee">'
                views = []
                for m in substr['matches']:
                    viewid += 1
                    viewerid = 'view'+str(viewid)
                    r1 = '<td>{:9f} in {}</br>{}</td>'.format(m['mz'], m['scan'], m['spectrum'])
                    r2 = '<td><script>show_mol("{}","{}","{}","{}")</script></td>'.format(viewerid, m['mol'].replace('\n','\\n'), m['fragatoms'], feature['type'])
                    views.append([r1, r2])
                for i in range(len(views)/6+1):
                    print '<tr>' + ''.join([m[0] for m in views[i*6:i*6+6]]) + '</tr>'
                    print '<tr>' + ''.join([m[1] for m in views[i*6:i*6+6]]) + '</tr>'
                print '</table></div></div></br>'
        print '<br>'

def write_collapsible_script():
    print """<script>
var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
}
</script>"""

# MAIN
motifannotations = json.load(sys.stdin)
write_header()
write_motif_annotations()
write_collapsible_script()

