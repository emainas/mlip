#!/usr/bin/env python3
"""Apply compact landscape/Avenir styling to a Pandoc-generated DOCX."""

from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile
import re
import shutil
import sys


def replace_fonts(xml: str) -> str:
    # Set both document defaults and named styles to Avenir while preserving
    # complex-script fallbacks declared by Word/Pandoc.
    xml = re.sub(
        r'<w:rFonts\b[^>]*/>',
        '<w:rFonts w:ascii="Avenir" w:hAnsi="Avenir" w:eastAsia="Avenir" w:cs="Avenir"/>',
        xml,
    )
    return xml


def style_document(xml: str) -> str:
    # US Letter landscape, 0.60-inch margins.
    xml = re.sub(
        r'<w:pgSz\b[^>]*/>',
        '<w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/>',
        xml,
    )
    xml = re.sub(
        r'<w:pgMar\b[^>]*/>',
        '<w:pgMar w:top="864" w:right="864" w:bottom="864" w:left="864" '
        'w:header="432" w:footer="432" w:gutter="0"/>',
        xml,
    )
    return xml


def main(src: Path, dst: Path) -> None:
    with TemporaryDirectory() as td:
        root = Path(td)
        with ZipFile(src) as zf:
            zf.extractall(root)

        styles = root / "word" / "styles.xml"
        styles.write_text(replace_fonts(styles.read_text()), encoding="utf-8")

        doc = root / "word" / "document.xml"
        doc.write_text(style_document(doc.read_text()), encoding="utf-8")

        # Reinforce Avenir in theme font slots used by some Word versions.
        theme = root / "word" / "theme" / "theme1.xml"
        if theme.exists():
            text = theme.read_text()
            text = re.sub(r'typeface="(?:Calibri|Cambria|Aptos|Times New Roman)"', 'typeface="Avenir"', text)
            theme.write_text(text, encoding="utf-8")

        with ZipFile(dst, "w", ZIP_DEFLATED) as zf:
            for path in sorted(root.rglob("*")):
                if path.is_file():
                    zf.write(path, path.relative_to(root))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: style_osteocalcin_docx.py INPUT.docx OUTPUT.docx")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
