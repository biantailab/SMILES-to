import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def convert_smiles(smiles_string, output_file, output_format='pdb', generate_3d=True, image_size=200):
    """
    Universal SMILES conversion function
    Args:
        smiles_string: SMILES string to convert
        output_file: Output file path
        output_format: 'pdb', 'mol', or 'png'
        generate_3d: Whether to generate 3D coordinates (for pdb/mol)
        image_size: Size for PNG images (default: 200px square)
    """
    print(f"Converting SMILES to {output_format.upper()}: {smiles_string}")
    
    try:
        print(f"Parsing SMILES: {smiles_string}")
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print("Invalid SMILES string")
            return False
        
        print(f"Successfully parsed SMILES")
        print(f"Molecular formula: {Chem.MolToSmiles(mol)}")
        print(f"Number of atoms: {mol.GetNumAtoms()}")
        
        if output_format == 'png':
            # For image generation, use 2D coordinates
            print(f"Generating 2D molecular image...")
            AllChem.Compute2DCoords(mol)
            
            print(f"Saving as PNG image: {output_file} ({image_size}x{image_size}px)")
            img = Draw.MolToImage(mol, size=(image_size, image_size))
            img.save(output_file)
        else:
            # For PDB/MOL files, add hydrogens and handle 3D
            mol_with_h = Chem.AddHs(mol)
            print(f"Number of atoms after adding hydrogens: {mol_with_h.GetNumAtoms()}")
            
            if generate_3d:
                print("Generating 3D coordinates...")
                try:
                    embed_result = AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
                    if embed_result == -1:
                        print("First embedding failed, trying alternative method...")
                        AllChem.EmbedMolecule(mol_with_h)
                    
                    minimize_result = AllChem.MMFFOptimizeMolecule(mol_with_h)
                    if minimize_result != 0:
                        print("Energy minimization may not have fully converged")
                    
                    print("3D coordinates generated successfully")
                except Exception as e:
                    print(f"3D coordinate generation warning: {e}")
                    print("Using 2D coordinates...")
                    AllChem.Compute2DCoords(mol_with_h)
            else:
                AllChem.Compute2DCoords(mol_with_h)
            
            # Save file based on format
            if output_format == 'pdb':
                print(f"Saving as PDB file: {output_file}")
                Chem.MolToPDBFile(mol_with_h, output_file)
            elif output_format == 'mol':
                print(f"Saving as MOL file: {output_file}")
                Chem.MolToMolFile(mol_with_h, output_file)
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            print(f"Conversion successful: {output_file}")
            return True
        else:
            print("Output file creation failed")
            return False
            
    except Exception as e:
        print(f"Conversion failed: {e}")
        return False

def interactive_conversion():
    """Interactive molecular conversion tool"""
    print("=" * 60)
    print("🔬 Interactive Molecular Format Conversion Tool")
    print("=" * 60)
    print("This tool converts SMILES strings to PDB, MOL, or PNG files")
    print("=" * 60)
    
    while True:
        print("\n" + "-" * 60)
        print("Options:")
        print("1. Convert SMILES to PDB (3D)")
        print("2. Convert SMILES to PDB (2D)")
        print("3. Convert SMILES to MOL (3D)")
        print("4. Convert SMILES to MOL (2D)")
        print("5. Convert SMILES to PNG (2D)")
        print("6. Exit")
        print("-" * 60)
        
        try:
            choice = input("Please select an option (1-6): ").strip()
            
            if choice == '6':
                print("\n👋 Goodbye!")
                break
            elif choice not in ['1', '2', '3', '4', '5']:
                print("❌ Invalid choice. Please select 1-6.")
                continue
            
            # Get SMILES input
            print("\n" + "=" * 40)
            smiles = input("Enter SMILES string: ").strip()
            if not smiles:
                print("❌ Please enter a valid SMILES string")
                continue
            
            # Determine format and settings based on choice
            if choice == '5':
                output_format = 'png'
                default_ext = '.png'
                generate_3d = False  # Images are always 2D
                
                # Get image size for PNG
                print("\n" + "-" * 40)
                print("Enter image size in pixels (e.g., 400 for 400x400px)")
                print("-" * 40)
                
                size_input = input("Image size (px) [default: 200px]: ").strip()
                if not size_input:
                    image_size = 200
                else:
                    try:
                        image_size = int(size_input)
                        if image_size <= 0:
                            print("❌ Invalid size, using default 200px")
                            image_size = 200
                    except ValueError:
                        print("❌ Invalid size, using default 200px")
                        image_size = 200
                    
            elif choice in ['1', '2']:
                output_format = 'pdb'
                default_ext = '.pdb'
                generate_3d = (choice == '1')
                image_size = 200  # Not used for PDB/MOL
            else:  # choice in ['3', '4']
                output_format = 'mol'
                default_ext = '.mol'
                generate_3d = (choice == '3')
                image_size = 200  # Not used for PDB/MOL
            
            # Get output filename
            output_name = input(f"Enter output filename (default: output{default_ext}): ").strip()
            if not output_name:
                output_name = f"output{default_ext}"
            elif not output_name.endswith(default_ext):
                output_name += default_ext
            
            print(f"\n🔄 Converting SMILES to {output_format.upper()}...")
            print("=" * 40)
            
            # Perform conversion using unified function
            success = convert_smiles(smiles, output_name, output_format, generate_3d, image_size)
            
            if success:
                print("\n✅ Operation completed successfully!")
                print(f"📁 Output file: {output_name}")
            else:
                print("\n❌ Operation failed!")
            
            # Ask if user wants to continue
            print("\n" + "-" * 40)
            continue_choice = input("Perform another operation? (y/n): ").strip().lower()
            if continue_choice not in ['y', 'yes']:
                print("\n👋 Goodbye!")
                break
                
        except KeyboardInterrupt:
            print("\n\n👋 Program interrupted. Goodbye!")
            break
        except Exception as e:
            print(f"\n❌ An error occurred: {e}")
            continue

def main():
    """Main function"""
    try:
        interactive_conversion()
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
