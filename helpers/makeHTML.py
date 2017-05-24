import pandas as pd;

def generateHTML(reportLoc,repoLoc,namesLoc,outputLoc):
    """
    generateHTML: generate HTML report from report.taxon: 
    
        reportLoc - location of report.taxon
        repoLoc - location of Pandora repository
        namesLoc - location of names.dmp file containing taxid to name
        outputLoc - location of file output folder
    """
    
    
    #reportLoc = "report.taxon.txt"
    canvasLoc = repoLoc+"/scripts/vendor/jquery.canvasjs.min.js"
    jqueryLoc = repoLoc+"/scripts/vendor/jquery-3.1.1.slim.min.js"
    datatablesLoc = repoLoc+"/scripts/vendor/datatables.min.js"
    #namesLoc = "/names.dmp"
    
    #load taxon report
    df = pd.read_table(reportLoc)
    df.applymap(str)
    df.columns = ['name','taxID','#_of_read','#_of_contigs','longest_contig','len_longest_contig','RPMH']
    sampleID = str(df['name'][0])
    
    #load canvasJS
    with open(canvasLoc,'r') as canvasjs:
        canvasString = canvasjs.read()
    
    #load jquery
    with open(jqueryLoc,'r') as jqueryjs:
        jqueryString = jqueryjs.read()
    
    #load DataTables
    with open(datatablesLoc,'r') as datatables:
        datatablesString = datatables.read()
    
    #load taxid to name 
    nameDump = pd.read_table(namesLoc,sep="|",header=None)
    nameDump[1] = nameDump[1].str.strip()
    nameDump[2] = nameDump[2].str.strip()
    nameDump[3] = nameDump[3].str.strip()
    nameDump = nameDump.applymap(str)
    #generate canvasJS string for donut chart
    pd.options.mode.chained_assignment = None
    donutString = ""
    for i in range(0,len(df.index)):
        df['taxID'][i] = df['taxID'][i].split(";")[0]
        df['name'][i] = nameDump[(nameDump[0] == df['taxID'][i]) & (nameDump[3] == 'scientific name')][[1]].to_string(header=False,index=False)
        temp = "{y:"+str(df['RPMH'][i])+",label:'"+df['name'][i]+" -'},"
        df['taxID'][i] = '<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id='+ str(df['taxID'][i]) +'">'+str(df['taxID'][i])+'</a>'
        donutString=donutString+temp
    donutString=donutString.rstrip(',')
    
    #workaround of column width to properly display taxid lnk
    old_width = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    tableString = df.to_html(index=False).replace("&lt;","<").replace("&gt;",">").replace('<table border="1" class="dataframe">','<table id="table1" border="1" class="dataframe">')
    pd.set_option('display.max_colwidth', old_width)
    
    #generate output including JS packages as string in header
    htmloutput = '''
    <html>
        <head >
            <style>
                h1,h2,table{
                    font-family: "Open Sans","Arial",sans-serif;
                    color:#3C3C3C;
                }
                h1{
                    font-size: 40px;
                }
                h2{
                    font-size: 28px;
                    padding-left: 5px;
                    margin-bottom: 5px;
                }
                body{
                    margin-left: 50px;
                    margin-right:50px;
                }
                tr:hover{
                    background-color:#f5f5f5;
                }
                th{
                    text-align:left;
                }
                table {
                    border-collapse:collapse;
                    border: none;
                    border-color:white;
                    border-bottom:1px solid #ddd;
                }
                .tablediv{
                    padding-left:10px;
                }
            </style>
            <link href="https://fonts.googleapis.com/css?family=Open+Sans" rel="stylesheet">
            <script type="text/javascript">'''\
            + jqueryString + '''
            </script>
            <script type="text/javascript">'''\
            + canvasString + '''
            </script>
            <script type="text/javascript">'''\
            + datatablesString + '''
            </script>
            <script type="text/javascript">
                $(document).ready( function () {
                    $('#table1').DataTable(
                    {
                        paging:false,
                        searching:false,
                        fixedHeader:true,
                        "order":[[6,'desc']],
                        scrollY: 300,
                        "info":false
                    }
                    );
                } );
            </script>
            <script type="text/javascript">
                window.onload = function () {
                    var chart = new CanvasJS.Chart(
                        "chartContainer",
                        {
                                    animationEnabled: true,
                            data: [
                            {        
                                type: "doughnut",
                                startAngle:20,
                                innerRadius: "60%",
                                toolTipContent: "{label}: {y} - <strong>#percent%</strong>",
                                indexLabel: "{label} #percent%",
                                dataPoints: [''' \
                                + donutString +'''
                                ]
                            }
                            ]
                        }
                    );
                    chart.render();
                }
            </script>
        </head>
        <body>
            <h1>PANDORA Scan - "''' + sampleID +'''"</h1>
            <h2>Pathogen Abundance</h2>
            <div id="chartContainer" style="height: 400px; width: 100%;"></div>
            <h2>Results</h2>
            <div class="tablediv">
            '''+ tableString +'''
            </div>
        </body>
    </html>'''
    
    #write to html file
    with open(outputLoc+'/report.taxon.html', 'w') as f:
        f.write(htmloutput)

