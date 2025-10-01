IMPORTANT: flag files need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. Collect all designs that passed the filtering steps in both 02_fastdesign and 03_grafting and split them into silent file batches using the silentsplitshuf command
2. Run the surface design RosettaScript ./surface_redesign.xml on each batch using ./fastdesign.flags
3. Collect all resulting designs into a single silent file