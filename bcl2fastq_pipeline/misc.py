"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
import xml.etree.ElementTree as ET
from reportlab.lib import colors
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from time import strftime
from reportlab.pdfgen import canvas

def makeProjectPDF(node, project, config) :
    """
    Uses reportlab to generate a per-project PDF file.

    Contents are currently a table containing:
      * Sample
      * Barcode
      * Lane
      * # Reads (passing filter)
      * % Bases > Q30
      * Average base quality
    For paired-end datasets, the last two columns are repeated and named differently.
    """
    stylesheet=getSampleStyleSheet()
    #stylesheet.add(ParagraphStyle(name="fixed", parent=stylesheet['Normal'], fontName="Courier", spaceBefore=10))

    pdf = BaseDocTemplate("%s/%s/%s/SequencingReport.pdf" % (
        config.get("Paths","outputDir"),config.get("Options","runID"),
        project), pagesize=A4)
    topHeight=100 #The image is 86 pixels tall
    fTL = Frame(pdf.leftMargin, pdf.height, width=pdf.width/2, height=topHeight, id="col1") #Fixed height
    fTR = Frame(pdf.leftMargin+pdf.width/2, pdf.height, width=pdf.width/2, height=topHeight, id="col2")
    fB = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height-topHeight, id="bottom")
    fM = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height, id="main")
    
    elements = []
    PE = False
    if(len(node[0][0][0][0][1]) == 3) :
        PE = True
        data = [["Sample ID","Barcode","Lane","# Reads","% Bases\n>= Q30\nRead #1","Ave. Qual.\nRead #1","% Bases\n>= Q30\nRead #2","Ave. Qual.\nRead #2"]]
    else :
        data = [["Sample ID","Barcode","Lane","# Reads","% Bases\n>= Q30","Ave. Qual."]]

    #A text blurb
    string = "Project: %s" % project
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "Report generated: %s" % (strftime("%d-%m-%Y %H:%M:%S"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "BCL2Fastq pipeline version: %s" % (config.get("Version","pipeline"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "bcl2fastq version: %s" % (config.get("Version","bcl2fastq"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    readLength = int(int(node[0][0][0][0][0][1][0].text)/int(node[0][0][0][0][0][0].text))
    string = "FastQC version: %s" % (config.get("Version","fastQC"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    if(PE) :
        string = "%i base paired-end reads" % readLength
    else :
        string = "%i base single-end reads" % readLength
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)

    #Image
    elements.append(Image(config.get("Options","imagePath"), hAlign="RIGHT"))
    elements.append(NextPageTemplate("RemainingPages"))

    #Data table
    for sample in node.findall("Sample") :
        e = [None,None,-1,0,0,0,0,0,0]
        e[0] = sample.get("name")
        if(e[0] == "all") :
            continue
        for barcode in sample.findall("Barcode") :
            e[1] = barcode.get("name")
            if(e[1] == "all") :
                continue
            for lane in barcode.findall("Lane") :
                e[2] = lane.get("number")
                e[3] = 0
                e[4] = 0
                e[5] = 0
                e[6] = 0
                e[7] = 0
                e[8] = 0
                for tile in lane.findall("Tile") :
                    e[3] += int(tile[1][0].text) #Pf->ClusterCount
                    e[4] += int(tile[1][1][0].text) #Pf->Read1->Yield
                    e[5] += int(tile[1][1][1].text) #Pf->Read1->YieldQ30
                    e[6] += int(tile[1][1][2].text) #Pf->Read1->QualSum
                    if(PE) :
                        e[7] += int(tile[1][2][1].text) #Pf->Read2->YieldQ30
                        e[8] += int(tile[1][2][2].text) #Pf->Read2->QualSum
                if(PE) :
                    data.append([e[0],
                                 e[1],
                                 e[2],
                                 e[3],
                                 "%5.2f" % (100*(e[5]/e[4])),
                                 "%5.2f" % (e[6]/e[4]),
                                 "%5.2f" % (100*(e[7]/e[4])),
                                 "%5.2f" % (e[8]/e[4])
                        ])
                else :
                    data.append([e[0],
                                 e[1],
                                 e[2],
                                 e[3],
                                 "%5.2f" % (100*(e[5]/e[4])),
                                 "%5.2f" % (e[6]/e[4])
                        ])

    t = Table(data, style=[
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), (0xD3D3D3, None)) #Light grey
        ], repeatRows=1)
    elements.append(t)

    #Add the key
    elements.append(Spacer(0,30))
    key = []
    key.append([Paragraph("Sample ID",
            stylesheet['BodyText']),
        Paragraph("The sample ID as provided to the sequencer in the sample sheet. This may not match the final file name, but will match the directory in which it's held.", 
            stylesheet['BodyText'])])
    key.append([Paragraph("Barcode",
            stylesheet['BodyText']),
        Paragraph("The sample barcode added by the sequencing facility (or you, if you created the libraries yourself). This will generally be 6 nucleotides long.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Lane", 
            stylesheet['BodyText']),
        Paragraph("The lane number on the flow cell (there are 8 of them).",
            stylesheet['BodyText'])])
    key.append([Paragraph("# Reads", 
            stylesheet['BodyText']),
        Paragraph("The number of reads in a given file. For paired-end datasets, this is equivalent to the number of fragments sequenced, rather than summing the counts for read #1 and read #2. Note that this includes only reads passing the quality filter.",
            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #1", 
            stylesheet['BodyText']),
        Paragraph("The percentage of bases in read #1 of a pair having a Phred-scaled score of at least 30, meaning that the 0.1% or less chance that they're incorrect.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual. Read #1", 
            stylesheet['BodyText']),
        Paragraph("The average Phred-scaled base quality of bases in read #1 of a pair. This number of -10*log10(Probability that the call is incorrect). In other words, if a call is 100% likely to be wrong, the score is 0 (or 10 for 10% likelihood, 20 for 1% likelihood, etc.).",
            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #2", 
            stylesheet['BodyText']),
        Paragraph("Identical to '% Bases >= Q30 Read #1', but for read #2 of a pair.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual. Read #2", 
            stylesheet['BodyText']),
        Paragraph("Identical to 'Ave. Qual. Read #1', but for read #1 of a pair.",
            stylesheet['BodyText'])])
    key.append([Paragraph("# Reads", 
            stylesheet['BodyText']),
        Paragraph("Identical to '% Bases >= Q30 Read #1', but for single-end datasets.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual.", 
            stylesheet['BodyText']),
        Paragraph("Identical to 'Ave. Qual. Read #1', but for single-end datasets.",
            stylesheet['BodyText'])])
    t2 = Table(key, colWidths=(80, None))
    t2.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP')]))
    elements.append(t2)

    pdf.addPageTemplates([PageTemplate(id="FirstPage", frames=[fTL, fTR, fB]),
        PageTemplate(id="RemainingPages", frames=[fM])]),
    pdf.build(elements)

def getFCmetrics(root) :
    barcode = root[0][0] #Sample "all", barcode "all"
    message = "Lane\t# Clusters (% pass)\t% Bases >=Q30\tAve. base qual.\n"
    for lane in barcode.findall("Lane") :
        message += "Lane %s" % lane.get("number")
        clusterCount = 0
        clusterCountPass = 0
        baseYield = [0,0]
        baseYieldQ30 = [0,0]
        QualSum = [0,0]
        rlens=[0,0]
        for tile in lane :
            clusterCount += int(tile[0][0].text)
            clusterCountPass += int(tile[1][0].text)
            #Yield
            baseYield[0] += int(tile[1][1][0].text)
            if(len(tile[1]) == 3) :
                baseYield[1] += int(tile[1][2][0].text)
            #YieldQ30
            baseYieldQ30[0] += int(tile[1][1][1].text)
            if(len(tile[1]) == 3) :
                baseYieldQ30[1] += int(tile[1][2][1].text)
            #QualSum
            QualSum[0] += int(tile[1][1][2].text)
            if(len(tile[1]) == 3) :
                QualSum[1] += int(tile[1][2][2].text)
        #Number of clusters (%passing filter)
        message += "\t%i (%5.2f%%)" % (clusterCount,100*clusterCountPass/clusterCount)
        #%bases above Q30
        if(baseYield[1] > 0) :
            message += "\t%5.2f%%/%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]),
                100*(baseYieldQ30[1]/baseYield[1]))
        else :
            message += "\t%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]))
        #Average base quality
        if(baseYield[1] > 0) :
            message += "\t%4.1f/%4.1f\n" % (QualSum[0]/clusterCountPass/(baseYield[0]/clusterCount),
                QualSum[1]/clusterCountPass/(baseYield[1]/clusterCount))
        else :
            message += "\t%4.1f\n" % (QualSum[0]/clusterCountPass/(baseYield[0]/clusterCount))

    return message

def parseConversionStats(config) :
    """
    Parse ConversionStats.xml, producing:
     1) A PDF file for each project
     2) A message that will be included in the email message
    """
    tree = ET.parse("%s/%s/Stats/ConversionStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    root = tree.getroot()[0] #We only ever have a single flow cell
    #Per-project PDF files
    for project in root.findall("Project") :
        if(project.get("name") == "default") :
            continue
        if(project.get("name") == "all") :
            metrics = getFCmetrics(project)
        else :
            makeProjectPDF(project, project.get("name"), config)
    return metrics

def enoughFreeSpace(config) :
    """
    Ensure that outputDir has at least minSpace gigs
    """
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    free /= 1024*1024*1024
    if(free >= float(config.get("Options","minSpace"))) :
        return True
    return False

def errorEmail(config, msg) :
    msg = MIMEText(msg)
    msg['Subject'] = "[bcl2fastq_pipeline] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()

def finishedEmail(config, msg, runTime) :
    message = "Flow cell: %s\n" % config.get("Options","runID")
    message += "Run time: %s\n" % runTime
    message += msg

    msg = MIMEText(message)
    msg['Subject'] = "[bcl2fastq_pipeline] %s processed" % config.get("Options","runID")
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
