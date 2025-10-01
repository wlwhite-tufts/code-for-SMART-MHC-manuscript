IMPORTANT: flag files need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files. Both RIFdock flag files make use of the xform_pos_ang10_0.35A_1.1d.x in this directory.

Part 1: RIFgen
	1. Select target residues
		a. use pymol to select MHC residues with CA-to-CA distance <= 9A to B2m (using the unrelaxed 1S7U structure), minus any pointing away from the interface
		b. write these residue numbers in ./rifgen_residues.list
	2. Generate standard RIF
		a. use standard rifgen flags (./standard_rifgen.flag), modified to reference the target (../00_prep/1S7U_relaxed_head.pdb) and the rifgen residues (./rifgen_residues.list)
		b. run rifgen (https://github.com/rifdock/rifdock/tree/master) using the specified flags (./standard_rifgen.flag)
	3. Select native B2m residues to generate "native" RIF
		a. use pymol to evaluate B2m residues in significant contact with a1+a2 domains
		b. delete all other residues from unrelaxed 1S7U structure and save as ./native_hotspots.pdb
	4. Generate the "native" RIF
		a. use modified flags (./native_rifgen.flag) to generate a RIF from ./native_hotspots.pdb
		b. target and RIFgen software are the same as the standard RIFgen step
	Note: RIF residues, and native hotspots are labeled in ./1S7U_annot.pse as rifres and nat_res, respectively

Part 2: PatchDock

	1. Follow steps 6 - 8 in the Cao et al guide
		a. prepare the target (center and relabel to chain B) 
		b. prepare the scaffold library (convert to polyVal) - scaffold library can be found here:
			http://files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/scaffolds.tar.gz
		c. generate PatchDock commands (with the -target_res flag set to ./patchdock_residues.list)
			This list of residues is nearly identical to the rifgen residues, but leaves off some more preipheral target residues
	2. Run PatchDock commands (creates a set of starting positions for RIFdock to use)

Part 3: RIFdock

	1. Setup standard RIF flags. Follow step 9 in the Cao et al guide to fill in the required information to the ./standard_rifdock.flag file using the rifgen.log file from part1:step2
	2. Setup native RIF flags. Repeat the setup for ./native_rifdock.flags, using the rifgen.log file from part1:step4
	3. Run rifdock on the standard RIF, following step 10 of the Cao et al guide, using ./standard_rifdock.flag provided here instead of tuning parameters as suggested by the guide.
	4. Run rifdock on the native RIF, following step 10 of the Cao et al guide, using ./native_rifdock.flag provided here instead of tuning parameters as suggested by the guide.