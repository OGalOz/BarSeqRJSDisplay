

How the JavaScript works:
Within the file brsq_viz.html, there is a section
where other Javascript files are imported. The basic
visualization files are all but the ones that 
are under 'VARIABLE Data'. Under 'VARIABLE Data',
there are 3 different files, each essentially a 
JSON object.
The first is called 'volcano_obj.js', which holds
a variable called 'volcano_object'. This object
has data that has the following keys:
            "min_y": float 
            "max_y": float
            "min_x": float
            "max_x": float
            "t_mega_matrix": list<num_list>
            "fit_mega_matrix": list<num_list>
where each 'num_list' is a list of floats - 
these lists simply correspond to the actual
matrices where the data begins.
The second is called 'exp2col', which has the 
following keys:
            "col2ix": dict, experiment + set name - > int representing column
            "ix2col": the reverse of above, but str(int) -> experiment + set name
The third is called 'genes2row', which has the 
following keys:
    "rownum2gene": str(int) -> gene info split by "|"
    "gene_id4_2rownum": gene info split by "|" w/o the last piece -> tuple<row_num (int), description (str)>
    "desc2id4": gene description (str) -> gene info split by "|" w/o the last piece
