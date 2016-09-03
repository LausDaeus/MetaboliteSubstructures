from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageChops


def draw_smiles(smiles_string, file_path, size=(1000,1000)):
    try:
        molecule = Chem.MolFromSmiles(smiles_string)
        Draw.MolToFile(molecule, file_path, size=size)
        trim_image(file_path)
    except ValueError:
        print "ERROR: Molecule could not be drawn to file."

def trim_image(file_path):
    im = Image.open(file_path)
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        im = im.crop(bbox)
        im.save(file_path, "PNG")

