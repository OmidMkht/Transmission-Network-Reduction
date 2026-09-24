"""Write CSVs into one multi-sheet .xlsx, using only the standard library.

    python analysis/xlsx_write.py out.xlsx SheetName=file.csv [SheetName=file.csv ...]

xlsx is a zip of XML. Strings go in inline (no shared-string table), which keeps
this to one function and is perfectly valid OOXML. Numeric-looking cells are
written as numbers so Excel can chart them; everything else is text. The first
row of each CSV becomes a bold header with a frozen pane and an autofilter.
"""
import csv
import os
import re
import sys
import zipfile
from xml.sax.saxutils import escape

NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def col_name(i):
    """1 -> A, 27 -> AA"""
    s = ""
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def sheet_xml(rows):
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">']
    ncol = max((len(r) for r in rows), default=1)
    out.append('<sheetViews><sheetView workbookViewId="0">'
               '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
               '</sheetView></sheetViews>')
    out.append('<sheetFormatPr defaultRowHeight="15"/>')
    out.append("<sheetData>")
    for ri, row in enumerate(rows, start=1):
        out.append(f'<row r="{ri}">')
        for ci, val in enumerate(row, start=1):
            ref = f"{col_name(ci)}{ri}"
            text = "" if val is None else str(val).strip()
            if text == "":
                continue
            if ri > 1 and NUM.match(text):
                out.append(f'<c r="{ref}"><v>{text}</v></c>')
            else:
                style = ' s="1"' if ri == 1 else ""
                out.append(f'<c r="{ref}"{style} t="inlineStr">'
                           f"<is><t>{escape(text)}</t></is></c>")
        out.append("</row>")
    out.append("</sheetData>")
    if rows:
        out.append(f'<autoFilter ref="A1:{col_name(ncol)}{len(rows)}"/>')
    out.append("</worksheet>")
    return "".join(out)


STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>
<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>
<fills count="2"><fill><patternFill patternType="none"/></fill>
<fill><patternFill patternType="gray125"/></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>
</styleSheet>"""


def write_xlsx(path, sheets):
    """sheets: list of (name, rows) where rows is a list of lists."""
    n = len(sheets)
    ct = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
          '<Default Extension="xml" ContentType="application/xml"/>',
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
    for i in range(1, n + 1):
        ct.append(f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
                  'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    ct.append("</Types>")

    wb = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
    for i, (name, _) in enumerate(sheets, start=1):
        safe = escape(name[:31].replace("/", "-").replace("\\", "-"))
        wb.append(f'<sheet name="{safe}" sheetId="{i}" r:id="rId{i}"/>')
    wb.append("</sheets></workbook>")

    rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for i in range(1, n + 1):
        rels.append(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/'
                    f'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>')
    rels.append(f'<Relationship Id="rId{n+1}" Type="http://schemas.openxmlformats.org/'
                'officeDocument/2006/relationships/styles" Target="styles.xml"/>')
    rels.append("</Relationships>")

    root = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/'
            '2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "".join(ct))
        z.writestr("_rels/.rels", root)
        z.writestr("xl/workbook.xml", "".join(wb))
        z.writestr("xl/_rels/workbook.xml.rels", "".join(rels))
        z.writestr("xl/styles.xml", STYLES)
        for i, (_, rows) in enumerate(sheets, start=1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", sheet_xml(rows))
    return path


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [r for r in csv.reader(fh)]


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    out = sys.argv[1]
    sheets = []
    for arg in sys.argv[2:]:
        name, _, path = arg.partition("=")
        if not path:
            name, path = os.path.splitext(os.path.basename(arg))[0], arg
        if not os.path.isfile(path):
            print(f"  skip {name}: {path} not found")
            continue
        sheets.append((name, read_csv(path)))
    if not sheets:
        print("nothing to write")
        sys.exit(1)
    write_xlsx(out, sheets)
    print(f"wrote {out}  ({len(sheets)} sheets: " +
          ", ".join(n for n, _ in sheets) + ")")
