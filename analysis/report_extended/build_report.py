from pathlib import Path
import json, re
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TCGA_Drug_Analysis_EN_20261013'
content=json.loads((ROOT/'analysis/report_extended/report_content.json').read_text())
assert content.get('actual_results_ready') is True, 'Final verified analysis inputs required'
OUT.mkdir(parents=True,exist_ok=True)
doc=Document()
sec=doc.sections[0]
sec.page_width=Inches(8.5); sec.page_height=Inches(11)
sec.top_margin=Inches(.72); sec.bottom_margin=Inches(.68)
sec.left_margin=Inches(.75); sec.right_margin=Inches(.75)
sec.footer_distance=Inches(.24)
for name in ['Normal','Body Text','Title','Subtitle','Heading 1','Heading 2','Caption']:
    s=doc.styles[name]; s.font.name='Calibri'; s.font.color.rgb=RGBColor(0,0,0)
    s.element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'Calibri')
for x in doc.styles.element.xpath('.//w:pBdr'): x.getparent().remove(x)
normal=doc.styles['Normal'];normal.font.size=Pt(11)
normal.paragraph_format.space_after=Pt(6)
normal.paragraph_format.line_spacing=1.08
normal.paragraph_format.widow_control=True
doc.styles['Title'].font.size=Pt(20);doc.styles['Title'].font.bold=True
doc.styles['Title'].paragraph_format.space_after=Pt(9)
doc.styles['Subtitle'].font.size=Pt(10)
doc.styles['Subtitle'].paragraph_format.space_after=Pt(10)
for name,size in [('Heading 1',14),('Heading 2',12)]:
    s=doc.styles[name];s.font.size=Pt(size);s.font.bold=True
    s.paragraph_format.space_before=Pt(11);s.paragraph_format.space_after=Pt(6)
    s.paragraph_format.keep_with_next=True
doc.styles['Caption'].font.size=Pt(9)
doc.styles['Caption'].paragraph_format.line_spacing=1.03
doc.styles['Caption'].paragraph_format.space_after=Pt(8)

def hyperlink(p,label,url):
    h=OxmlElement('w:hyperlink')
    h.set(qn('r:id'),p.part.relate_to(url,RELATIONSHIP_TYPE.HYPERLINK,is_external=True))
    r=OxmlElement('w:r');pr=OxmlElement('w:rPr')
    color=OxmlElement('w:color');color.set(qn('w:val'),'1F497D');pr.append(color)
    u=OxmlElement('w:u');u.set(qn('w:val'),'single');pr.append(u)
    r.append(pr);t=OxmlElement('w:t');t.text=label;r.append(t);h.append(r);p._p.append(h)

def inline(p,text):
    # Source markdown links become native clickable Word hyperlinks.
    pat=re.compile(r'(\[([^\]]+)\]\((https?://(?:[^\s()]|\([^)]*\))+)\)|\*\*([^*]+)\*\*)')
    pos=0
    for m in pat.finditer(text):
        p.add_run(text[pos:m.start()])
        if m.group(2):hyperlink(p,m.group(2),m.group(3))
        else:p.add_run(m.group(4)).bold=True
        pos=m.end()
    p.add_run(text[pos:])

def para(text,style=None):
    p=doc.add_paragraph(style=style);inline(p,text);return p

def table(block):
    headers=block['headers'];rows=block['rows'];widths=block.get('widths',[7/len(headers)]*len(headers))
    t=doc.add_table(rows=1,cols=len(headers));t.autofit=False;t.alignment=WD_TABLE_ALIGNMENT.CENTER
    borders=OxmlElement('w:tblBorders')
    for edge in ['top','left','bottom','right','insideH','insideV']:
        e=OxmlElement('w:'+edge);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'D9D9D9');borders.append(e)
    t._tbl.tblPr.append(borders)
    for i,w in enumerate(widths):t.columns[i].width=Inches(w)
    for i,h in enumerate(headers):t.rows[0].cells[i].text=str(h)
    for vals in rows:
        cells=t.add_row().cells
        for cell,val in zip(cells,vals):cell.text='' if val is None else str(val)
    for ri,row in enumerate(t.rows):
        trpr=row._tr.get_or_add_trPr();trpr.append(OxmlElement('w:cantSplit'))
        if ri==0:
            x=OxmlElement('w:tblHeader');x.set(qn('w:val'),'true');trpr.append(x)
        for ci,cell in enumerate(row.cells):
            cell.width=Inches(widths[ci]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            tcpr=cell._tc.get_or_add_tcPr();mar=OxmlElement('w:tcMar')
            for k,v in [('top',70),('bottom',70),('left',90),('right',90)]:
                e=OxmlElement('w:'+k);e.set(qn('w:w'),str(v));e.set(qn('w:type'),'dxa');mar.append(e)
            tcpr.append(mar)
            if ri==0:
                sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E2E2E2');tcpr.append(sh)
            for p in cell.paragraphs:
                p.paragraph_format.space_after=Pt(0);p.paragraph_format.line_spacing=1.02
                if ci in block.get('center_columns',[]):p.alignment=WD_ALIGN_PARAGRAPH.CENTER
                for r in p.runs:r.font.size=Pt(block.get('font_size',9.5));r.bold=(ri==0)
    doc.add_paragraph().paragraph_format.space_after=Pt(0)

def figure(block):
    path=ROOT/block['path']
    if not path.exists():raise FileNotFoundError(path)
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next=True
    p.add_run().add_picture(str(path),width=Inches(block.get('width',6.7)))
    # Accessibility descriptions are separate from figure captions.
    for el in p._p.xpath('.//wp:docPr'):el.set('descr',block['caption'])
    para(block['caption'],'Caption')

footer=sec.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.RIGHT
fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');footer._p.append(fld)
doc.add_paragraph(content['title'],'Title')
doc.add_paragraph(content['subtitle'],'Subtitle')
for b in content['blocks']:
    kind=b['type']
    if kind=='paragraph':para(b['text'])
    elif kind=='heading':para(re.sub(r'[^A-Za-z0-9 ]',' ',b['text']),'Heading '+str(b.get('level',1)))
    elif kind=='table':table(b)
    elif kind=='figure':figure(b)
    elif kind=='pagebreak':doc.add_page_break()
    elif kind=='reference':
        p=para(b['text']);p.paragraph_format.space_after=Pt(4);p.paragraph_format.line_spacing=1.0
        for r in p.runs:r.font.size=Pt(9.5)
    else:raise ValueError(kind)
doc.core_properties.title=content['title']
doc.core_properties.subject='TCGA and drug resource analysis of CRISPR olaparib resistance candidates'
doc.core_properties.author=''
doc.save(OUT/'Report_EN.docx')
print(OUT/'Report_EN.docx')
