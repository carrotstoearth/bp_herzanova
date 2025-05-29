# Gene Dotmap Pipeline

Tato pipeline slouží k analýze a vizualizaci výskytu genů indikujících horizontální (HGT) a vertikální (VGT) genový přenos napříč souborem bakteriálních genomů. Pomocí tří propojených Python skriptů lze:

- zpracovat GFF anotace,
- vytvořit porovnávací matice přítomnosti genů,
- vizualizovat sdílení genů ve formě statické i interaktivní dotmapy.

## 📁 Struktura

- `make_comparison_matrices.py` — vytvoří přítomnostní a porovnávací matice genů podle fylogenetického stromu.
- `generate_dotmaps.py` — vygeneruje kombinovanou statickou dotmapu všech zadaných genů.
- `dotmap_viewer.py` — spustí GUI aplikaci pro interaktivní vizualizaci výskytu HGT a VGT genů.

---

## 🧪 1. Generování porovnávacích matic

```bash
python make_comparison_matrices.py   --input_folder path/to/genome_folders   --output_folder path/to/output_dir   --tree_file path/to/tree_file.newick
```

### Vstupy:
- `--input_folder`: složka obsahující podadresáře s GFF/FFN soubory (např. výstup z Prokka).
- `--output_folder`: cílová složka pro výstupy.
- `--tree_file`: fylogenetický strom ve formátu Newick, podle kterého budou matice seřazeny.

### Výstupy:
- binární matice přítomnosti genů (CSV),
- porovnávací matice (`*_comparison_matrix_ordered_by_tree.csv`) pro každý HGT a VGT gen.

---

## 🖼️ 2. Generování statické dotmapy

```bash
python generate_dotmaps.py   --hgt_folder path/to/output_dir/hgt   --vgt_folder path/to/output_dir/vgt
```

### Popis:
Skript vykreslí kombinovanou dotmapu pro všechny geny ve složkách `hgt/` a `vgt/`. Každý bod odpovídá sdílení konkrétního genu mezi dvěma genomy.

### Výstupy:
- `combined_dotmap.png` uložený ve složce `hgt/` a `vgt/`.

---

## 🧭 3. Interaktivní vizualizace

```bash
python dotmap_viewer.py   --hgt_folder path/to/output_dir/hgt   --vgt_folder path/to/output_dir/vgt   --phylum_file path/to/phylum_metadata.xlsx
```

### Popis:
Spustí interaktivní GUI aplikaci pro výběr HGT a VGT genů a jejich vizualizaci ve formě dotmapy. Uživatel si může vybrat geny a zobrazit jejich výskyt napříč genomovou kolekcí.

### Parametry:
- `--hgt_folder`: složka obsahující CSV matice pro HGT geny.
- `--vgt_folder`: složka obsahující CSV matice pro VGT geny.
- `--phylum_file`: (volitelné) tabulka `.xlsx` s informací o fylu pro každý genom (sloupce „Strain“, „GTDB-Tk“).

---

## 🔧 Závislosti

Instalace požadovaných balíčků přes `pip`:

```bash
pip install biopython pandas numpy matplotlib PyQt5 plotly openpyxl
```

### Použité knihovny:
- `pandas`, `numpy`, `matplotlib`
- `PyQt5`, `plotly`, `openpyxl`
- `biopython` (pro práci s fylogenetickým stromem)

---

## 📂 Doporučená struktura vstupních dat

```
project_root/
├── genomes/
│   ├── genome1/
│   │   ├── genome1.gff
│   │   └── genome1.ffn
│   ├── genome2/
│   │   └── ...
├── tree_file.newick
├── phylum_metadata.xlsx
```

---

## 📌 Poznámky

- Seznam HGT a VGT genů je definován přímo ve skriptu `make_comparison_matrices.py` (proměnné `hgt_genes` a `vgt_genes`). Lze upravit podle potřeby.
- Skripty byly testovány na datech získaných nástroji Prokka a UBCG.
