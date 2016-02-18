- Some figures in PDF were converted to EPS with the following:

http://tex.stackexchange.com/questions/20883/how-to-convert-pdf-to-eps

Then edited BoundingBox by hand to ensure enough margin on the sides.  If new PDF plots are made:

1.)  SAVE THE OLDER EPS IN THIS FOLDER SOMEWHERE
2.)  Convert to EPS
3.)  Open EPS in text editor, update BoundingBox to match older EPS (should be sufficient, tweak if needed).

----------------------------------

You can compile with:
  latex gphoton.tex
  bibtex gphoton
  latex gphoton.tex

Assuming you have installed and upgraded the texlive package.

----------------------------------
