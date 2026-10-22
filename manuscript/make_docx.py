"""Build a Word (.docx) copy of the manuscript from the LaTeX sources for coauthors who edit in Word.

Usage (from manuscript/):  python3 make_docx.py [output.docx]

Pipeline: LaTeX sections -> pandoc markdown -> pandoc docx with citeproc (Springer Vancouver style, references.bib).
Equations become native Word equations, figure PDFs are rasterized to PNG, and Word line numbering is added.
Needs: pandoc, pdftoppm. The docx is a snapshot: re-run after editing the LaTeX (edits made in Word are not synced back).
"""
import re, subprocess, sys, tempfile, zipfile, shutil
from pathlib import Path

M = Path(__file__).resolve().parent
out = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else M.parent / "pycea_manuscript.docx"
MONTHS = ["January","February","March","April","May","June","July","August","September","October","November","December"]
rd = lambda p: (M / p).read_text()

def clean(s):
    s = re.sub(r"(?m)^%.*\n", "", s)
    s = re.sub(r"\\code\{([^{}]*)\}", lambda m: "\\texttt{" + m.group(1) + "}", s)
    return re.sub(r"\\label\{[^}]*\}", "", s)

tmp = Path(tempfile.mkdtemp(prefix="pycea_docx_"))
(tmp / "img").mkdir()
for n in ["1", "2", "3", "4", "S1"]:
    subprocess.run(["pdftoppm", "-r", "220", "-png", "-singlefile", str(M / f"figures/figure{n}.pdf"), str(tmp / f"img/figure{n}")], check=True)

# bibliography for citeproc: webpages -> online, access dates -> urldate, titles case-protected
bib = rd("references.bib").replace("@webpage{", "@online{")
bib = re.sub(r"lastchecked=\{(\d+) (\w+) (\d+)\}", lambda m: "urldate={%s-%02d-%02d}" % (m.group(3), MONTHS.index(m.group(2)) + 1, int(m.group(1))), bib)
bib = re.sub(r"(?m)^  title=\{(.*)\},$", lambda m: "  title={{" + m.group(1) + "}},", bib)
(tmp / "refs.bib").write_text(bib)

# title block from main.tex
mt = rd("main.tex")
title = re.search(r"\\title\{(.*)\}", mt).group(1)
kw = re.search(r"\\textbf\{Keywords:\}\s*(.*)", mt).group(1).strip()
authors = re.findall(r"\\author\[([^\]]*)\]\{([^}]*)\}", mt)
affils = re.findall(r"\\affil\[([^\]]*)\]\{(.*)\}", mt)
date = re.search(r"\\date\{\\small (.*)\}", mt).group(1)
sup = lambda t: "\\textsuperscript{%s}" % t
auth_line = ", ".join("\\textbf{%s}%s" % (n.replace("~", " "), sup(a)) for a, n in authors[:-1])
auth_line += ", and \\textbf{%s}%s" % (authors[-1][1].replace("~", " "), sup(authors[-1][0]))
aff_lines = " \\\\\n".join(sup(("*" if k == "*" else k)) + v for k, v in affils)
head = auth_line + "\n\n" + aff_lines + "\n\n" + date + "\n"

fig = clean(rd("sections/figures.tex"))
fig = re.sub(r"\\includegraphics\[[^\]]*\]\{figure(\w+)\.pdf\}", lambda m: "\\includegraphics[width=\\textwidth]{img/figure%s.png}" % m.group(1), fig)
nums = iter(["1", "2", "3", "4", "S1"])
fig = re.sub(r"\\caption\{\\textbf\{", lambda m: "\\caption{\\textbf{Figure %s. " % next(nums), fig)
fig = re.sub(r"\\setcounter\{figure\}\{0\}\n", "", fig)
fig = re.sub(r"\\renewcommand\{\\thefigure\}\{[^\n]*\n", "", fig)
fig = fig.replace("{\\large\\textbf{Supplementary Information}}\\par\\medskip", "")
fig = re.sub(r"\\begin\{figure\}\[[a-z!]*\]", lambda m: "\\begin{figure}", fig)
i = fig.rindex("\\begin{figure}")
fig = fig[:i] + "\\section*{Supplementary Information}\n\n" + fig[i:]

body = "\\documentclass{article}\n\\begin{document}\n" + head + "\n\\section*{Abstract}\n" + clean(rd("sections/abstract.tex")) + "\n\n\\textbf{Keywords:} " + kw + "\n"
for t, f in [("Background", "introduction"), ("Results", "results"), ("Discussion", "discussion"), ("Conclusions", "conclusions"), ("Methods", "methods"), ("Declarations", "declarations")]:
    body += "\n\\section*{%s}\n" % t + clean(rd(f"sections/{f}.tex")) + "\n"
body += "\n\\section*{References}\n\nREFERENCES_PLACEHOLDER\n\n" + fig + "\n\\end{document}\n"
(tmp / "doc.tex").write_text(body)

md = subprocess.run(["pandoc", "-f", "latex", "-t", "markdown", "--wrap=none", "doc.tex"], cwd=tmp, check=True, capture_output=True, text=True).stdout
md = md.replace("REFERENCES_PLACEHOLDER", "::: {#refs}\n:::")
(tmp / "doc.md").write_text("---\ntitle: \"%s\"\n---\n\n" % title + md)
subprocess.run(["pandoc", "doc.md", "-f", "markdown+tex_math_dollars", "-o", "raw.docx", "--citeproc", "--bibliography=refs.bib",
                f"--csl={M / 'springer-vancouver-brackets.csl'}", "--resource-path=."], cwd=tmp, check=True)

# page setup (US Letter, 1 inch margins) and continuous line numbers, as in the LaTeX version; declare png type
with zipfile.ZipFile(tmp / "raw.docx") as zin, zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == "word/document.xml":
            xml = data.decode()
            # drop pandoc's bookmarks (Word draws them as [ ] brackets) and the internal links that pointed at them
            xml = re.sub(r"<w:bookmarkStart[^>]*/>|<w:bookmarkEnd[^>]*/>", "", xml)
            xml = re.sub(r'<w:hyperlink w:anchor="[^"]*"[^>]*>(.*?)</w:hyperlink>',
                         lambda m: m.group(1).replace('<w:rStyle w:val="Hyperlink" />', ""), xml, flags=re.S)
            data = xml.replace("<w:sectPr />", '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/><w:lnNumType w:countBy="1" w:restart="continuous"/></w:sectPr>').encode()
        if item.filename == "[Content_Types].xml" and b'Extension="png"' not in data:
            data = data.replace(b'<Default Extension="xml"', b'<Default Extension="png" ContentType="image/png" /><Default Extension="xml"', 1)
        zout.writestr(item, data)
shutil.rmtree(tmp)
print("wrote", out)
