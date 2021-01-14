#!python3

import os
import sys
import logging
import json
import shutil

"""
We convert output files from BarSeqR into JSON files
    so they'll be readable by a visualization program.

Generally, fitness and t score TSV files have the following headers:

    orgId
    locusId
    sysName
    geneName
    desc

and then a variable number of sets with different conditions
"""

def create_new_HTML_dir(HTML_out_dir, inp_fit, inp_t, org_name="Unknown",
                        col_start="5"):
    """ We make an HTML directory from which to run visualization
    
        We copy the default 'HTML' dir to the HTML_out_dir location
        Then we process the input fitness and t score files and create
        Javascript files in the HTML_out_dir/JS/directory
        Finally we updated the organisms name
    """

    if not os.path.exists(HTML_out_dir):
        #if os.access(HTML_out_dir, os.W_OK | os.X_OK):
        shutil.copytree("HTML", HTML_out_dir)
        #else:
        #raise Exception("You do not have write and execute access to " + \
        #           "output dir. Please choose a new destination.")
    else:
        raise Exception("Output directory already exists. Please choose new destination")

    add_files_to_js_dir(HTML_out_dir, 
                       inp_fit,
                       inp_t,
                       col_start)

    add_strings_to_brsq_viz_html(org_name, HTML_out_dir)

    os.system("open " + os.path.join(HTML_out_dir, "brsq_viz.html"))



def add_strings_to_brsq_viz_html(org_name, HTML_out_dir):
    """We add a title and info regarding organism name
    """
    html_FH = open(os.path.join(HTML_out_dir, "brsq_viz.html"),"r")

    new_f_str = ''

    c_line = html_FH.readline()
    c_ix = 0
    while c_line != "":
        c_ix += 1
        if c_ix == 5:
            # title line
            c_line = '    <title>' + org_name + ' BarSeqR Visualization</title>\n' 
        elif c_ix == 65:
            # organism name
            c_line = '    window.BarSeqROrganismName = "' + org_name + '"\n'
        new_f_str += c_line
        c_line = html_FH.readline()

    html_FH.close()
    html_FH = open(os.path.join(HTML_out_dir, "brsq_viz.html"),"w")
    html_FH.write(new_f_str)
    html_FH.close()



def add_files_to_js_dir(op_dir, 
                        inp_fit_score_fp,
                        inp_t_score_fp,
                        column_start_index):
    """
    All inputs string, col_start is string of int.
    """
    
    JS_dir = os.path.join(op_dir, "JS")
    column_start_index = int(column_start_index)

    if not os.path.isdir(JS_dir):
        raise Exception("JS dir not found in HTML dir? Must redownload.")


    # Checking existence of input TSV files
    for x in [inp_fit_score_fp, inp_t_score_fp]:
        if not os.path.isfile(x):
            raise Exception("Input file {} does not exist".format(
                            x))


    Fit_tsv_fh = open(inp_fit_score_fp, "r")
    T_tsv_fh = open(inp_t_score_fp, "r")


    # We get the headers in lists
    fit_headers = Fit_tsv_fh.readline().rstrip().split('\t')
    t_headers = T_tsv_fh.readline().rstrip().split('\t')

    # Below fit_column_headers refers to condition sets
    fit_column_headers = fit_headers[column_start_index:]
    t_column_headers = t_headers[column_start_index:]
    if len(fit_column_headers) != len(t_column_headers):
        raise Exception("fitness score file and t score file not " + \
                "compatible: t - {}, fit - {}".format(
                    inp_fit_score_fp, inp_t_score_fp))
    else:
        for i in range(len(fit_column_headers)):
            if fit_column_headers[i] != t_column_headers[i]:
                raise Exception("fitness score file and t score file not " + \
                        "compatible: t - {}, fit - {}".format(
                            inp_fit_score_fp, inp_t_score_fp))

    # We switch the 'desc' part of the headers to have condition before set #
    for i in range(len(t_column_headers)):
        desc = t_column_headers[i]
        desc_list = desc.split(' ')
        new_desc_list = desc_list[1:] + desc_list[0:1]
        new_desc = ' '.join(new_desc_list)
        t_column_headers[i] = new_desc


    # We get map from index number to column name
    ix_to_clm_nm = {i:t_column_headers[i] for i in range(len(t_column_headers))}
    # We get their reverse map as well
    clm_nm_to_ix = {v: k for k, v in ix_to_clm_nm.items()}


    fit_mega_matrix = []
    row_num_to_gene_info = {}

    # Now we take all the rows from the fitness TSV file 
    # and add it to the fitness mega_matrix
    # while simultaneously getting a map of row number
    # to gene info.
    c_row_num = 0
    c_line = Fit_tsv_fh.readline()
    while c_line != "":
        c_line = c_line.rstrip()
        c_list = c_line.split('\t')

        start_list= c_list[:column_start_index]
        row_num_to_gene_info[c_row_num] = "|".join(start_list)

        mm_list = c_list[column_start_index:]
        mm_list = [float(x) for x in mm_list]
        fit_mega_matrix.append(mm_list)
        
        c_row_num += 1
        c_line = Fit_tsv_fh.readline()


    t_mega_matrix = []
    c_line = T_tsv_fh.readline()
    while c_line != "":
        c_line = c_line.rstrip()
        c_list = c_line.split('\t')

        mm_list = c_list[column_start_index:]
        mm_list = [abs(float(x)) for x in mm_list]
        t_mega_matrix.append(mm_list)

        c_line = T_tsv_fh.readline()


    
    # Getting min and max fitness scores and t scores 
    # and populating fit to t score list

    max_fit = 0
    min_fit = 0
    for i in range(len(fit_mega_matrix)):
        c_row = fit_mega_matrix[i]
        for j in range(len(c_row)):
            c_item = float(c_row[j])
            if c_item > max_fit:
                max_fit = c_item
            elif c_item < min_fit:
                min_fit = c_item


    max_t = 0
    min_t = 0
    for i in range(len(t_mega_matrix)):
        c_row = t_mega_matrix[i]
        for j in range(len(c_row)):
            c_item = float(c_row[j])
            if c_item > max_t:
                max_t = c_item
            elif c_item < min_t:
                min_t = c_item
   
    volcano_dict = {
            "min_y": min_t,
            "max_y": max_t,
            "min_x": min_fit,
            "max_x": max_fit,
            "t_mega_matrix": t_mega_matrix,
            "fit_mega_matrix": fit_mega_matrix
    }
    
    
    with open(os.path.join(JS_dir,'volcano_obj.js'), "w") as g:
        g.write("window.volcano_object = " + \
                json.dumps(volcano_dict))
    logging.info("Wrote Output Volcano dict to "  + \
                os.path.join(JS_dir,'volcano_obj.js'))

    barseqr_column_info_d = {
        "col2ix" : clm_nm_to_ix,
        "ix2col": ix_to_clm_nm
    }
    
    logging.info("Writing exps to column number JSON at " + \
            os.path.join(JS_dir,'exp2col.js'))
    with open(os.path.join(JS_dir,'exp2col.js'), "w") as g:
        g.write("window.exp2col = " + \
                json.dumps(barseqr_column_info_d))

    barseqr_row_info_d = {
        "rownum2gene": row_num_to_gene_info,
        "gene2rownum": {v: k for k, v in row_num_to_gene_info.items()}
    }

    logging.info("Writing genes to row number JSON at " + \
                os.path.join(JS_dir,'genes2row.js'))
    with open(os.path.join(JS_dir,'genes2row.js'), "w") as g:
        g.write("window.gene2row = " + \
                json.dumps(barseqr_row_info_d))

    logging.info("Wrote all output JS files.")

    return None


    



def convert_tsvs_to_json(
                        inp_fit_score_fp,
                        inp_t_score_fp,
                        op_fit_json,
                        op_t_json,
                        op_volcano_dict_json,
                        experiments_to_column_number_json,
                        genes_to_row_number_json,
                        column_start_index=5):
    """


    Args:
        column_start_index: (int) Refers to the location
                            at which the indexes are IT
                            as opposed to preliminary info.
                            Default is 5
        all other inputs are str
    """

    # Checking existence of input TSV files
    for x in [inp_fit_score_fp, inp_t_score_fp]:
        if not os.path.isfile(x):
            raise Exception("Input file {} does not exist".format(
                            x))


    Fit_tsv_fh = open(inp_fit_score_fp, "r")
    T_tsv_fh = open(inp_t_score_fp, "r")


    # We get the headers in lists
    fit_headers = Fit_tsv_fh.readline().rstrip().split('\t')
    t_headers = T_tsv_fh.readline().rstrip().split('\t')

    # Below fit_column_headers refers to condition sets
    fit_column_headers = fit_headers[column_start_index:]
    t_column_headers = t_headers[column_start_index:]
    if len(fit_column_headers) != len(t_column_headers):
        raise Exception("fitness score file and t score file not " + \
                "compatible: t - {}, fit - {}".format(
                    inp_fit_score_fp, inp_t_score_fp))
    else:
        for i in range(len(fit_column_headers)):
            if fit_column_headers[i] != t_column_headers[i]:
                raise Exception("fitness score file and t score file not " + \
                        "compatible: t - {}, fit - {}".format(
                            inp_fit_score_fp, inp_t_score_fp))

    # We get map from index number to column name
    ix_to_clm_nm = {i:t_column_headers[i] for i in range(len(t_column_headers))}
    # We get ther reverse map as well
    clm_nm_to_ix = {v: k for k, v in ix_to_clm_nm.items()}


    fit_mega_matrix = []
    row_num_to_gene_info = {}

    # Now we take all the rows from the fitness TSV file 
    # and add it to the fitness mega_matrix
    # while simultaneously getting a map of row number
    # to gene info.
    c_row_num = 0
    c_line = Fit_tsv_fh.readline()
    while c_line != "":
        c_line = c_line.rstrip()
        c_list = c_line.split('\t')

        start_list= c_list[:column_start_index]
        row_num_to_gene_info[c_row_num] = "|".join(start_list)

        mm_list = c_list[column_start_index:]
        mm_list = [float(x) for x in mm_list]
        fit_mega_matrix.append(mm_list)
        
        c_row_num += 1
        c_line = Fit_tsv_fh.readline()


    t_mega_matrix = []
    c_line = T_tsv_fh.readline()
    while c_line != "":
        c_line = c_line.rstrip()
        c_list = c_line.split('\t')

        mm_list = c_list[column_start_index:]
        mm_list = [abs(float(x)) for x in mm_list]
        t_mega_matrix.append(mm_list)

        c_line = T_tsv_fh.readline()


    
    # Getting min and max fitness scores and t scores 
    # and populating fit to t score list

    max_fit = 0
    min_fit = 0
    for i in range(len(fit_mega_matrix)):
        c_row = fit_mega_matrix[i]
        for j in range(len(c_row)):
            c_item = float(c_row[j])
            if c_item > max_fit:
                max_fit = c_item
            elif c_item < min_fit:
                min_fit = c_item


    max_t = 0
    min_t = 0
    for i in range(len(t_mega_matrix)):
        c_row = t_mega_matrix[i]
        for j in range(len(c_row)):
            c_item = float(c_row[j])
            if c_item > max_t:
                max_t = c_item
            elif c_item < min_t:
                min_t = c_item
   
    volcano_dict = {
            "min_y": min_t,
            "max_y": max_t,
            "min_x": min_fit,
            "max_x": max_fit,
            "t_mega_matrix": t_mega_matrix,
            "fit_mega_matrix": fit_mega_matrix
    }
    
    
    with open(op_volcano_dict_json, "w") as g:
        g.write("window.volcano_object = " + \
                json.dumps(volcano_dict, indent=4))
    logging.info("Wrote Output Volcano dict to " + op_volcano_dict_json)

    barseqr_column_info_d = {
        "col2ix" : clm_nm_to_ix,
        "ix2col": ix_to_clm_nm
            }
    
    logging.info("Writing exps to column number JSON at " + \
            experiments_to_column_number_json)
    with open(experiments_to_column_number_json, "w") as g:
        g.write("window.exp2col = " + \
                json.dumps(barseqr_column_info_d, indent=4))

    barseqr_row_info_d = {
        "rownum2gene": row_num_to_gene_info,
        "gene2rownum": {v: k for k, v in row_num_to_gene_info.items()}
    }

    logging.info("Writing genes to row number JSON at " + \
                genes_to_row_number_json)
    with open(genes_to_row_number_json, "w") as g:
        g.write("window.gene2row = " + \
                json.dumps(barseqr_row_info_d, indent=4))


    



def convert_tsv_to_jsonOLD(inp_tsv_fp, op_json, 
                        score_map_json, sorted_scores_json):
    """ Take input TSV of fitness or t scores and convert to JSON

    Args:
        inp_tsv_fp: (str) path to TSV file
        op_json: (str) path to JSON file to write to
        score_map_json: (str) Path to JSON file which contains
            the score mapping.
            score -> list<loc1, loc2, ...>
    """

    if not os.path.isfile(inp_tsv_fp):
        raise Exception("Input file {} does not exist".format(
                        inp_tsv_fp))
    if os.path.isfile(op_json):
        raise Exception("Output file {} already exists. No overwrite.".format(
                        op_json))
    tsv_fh = open(inp_tsv_fp, "r")
    # We get the headers in a list
    headers = tsv_fh.readline().rstrip().split('\t')
    # We create a dict from header to loc
    header_d = {}
    for i in range(len(headers)):
        header_d[headers[i]] = i
        
    # Below list will contain all rows of file as lists
    tsv_list = []
    all_scores_set = set()

    c_line = tsv_fh.readline()
    while c_line != "":
        c_line = c_line.rstrip()
        c_list = c_line.split('\t')
        # converting float values to float
        for i in range(5,len(c_list)):
            c_list[i] = float(c_list[i])
            all_scores_set.add(c_list[i])
        tsv_list.append(c_list)
        c_line = tsv_fh.readline()

    tsv_fh.close()

    score_map = create_score_map_json(tsv_list)

    all_scores_list = sorted(list(all_scores_set))

    index_to_val = {header_d[s]:s for s in header_d.keys()}

    complete_json = {
            "value_to_index": header_d,
            "index_to_value": index_to_val,
            "list_of_values": tsv_list
    }

    with open(op_json, "w") as g:
        g.write(json.dumps(complete_json, indent=2))

    with open(score_map_json, "w") as g:
        g.write(json.dumps(score_map, indent=2))

    with open(sorted_scores_json, "w") as g:
        g.write(json.dumps(all_scores_list))

    logging.info("Wrote JSON files to {} and to {}".format(
        op_json, score_map_json))

    
def create_score_map_json(tsv_list):
    """We create an ordering of all the scores and map them to location

    In detail, we want to sort all the scores from lowest to highest,
    and for each score, note the location/s in which it is found, e.g.
        score 0.48 is found at [215,35] & [1246, 42]
        meaning that it is in the tsv_list at the 215th index list,
            and inside that list at the 35th index,
            AND
            at the tsv_list at the 1246th index, and the 42nd index
            of that list.
    Note that the beginning indeces of the values starts at the 6th 
        index (5 in 0-indexing) of every sublist of the input tsv_list.

    Args:
        tsv_list: list<list<str, str, str, str, str, float, ...>>
    Returns:
        score_map:
            score (float) -> list<location>
            where location: (list <row index, column index>)
            e.g. 
                {
                0.48: [[215, 35], [1246, 42]],
                .
                .
                .
                3.2: [[357, 42]]
                }
    """
    score_map = {}
    for k in range(len(tsv_list)):
        sublist = tsv_list[k]
        for j in range(5, len(sublist)):
            if sublist[j] in score_map:
                score_map[sublist[j]].append([k,j])
            else:
                score_map[sublist[j]] = [[k,j]]

    return score_map


def main():
    '''
    # Input takes HTML_dir output location,
    # an organism name, and the two TSV files
    # Then it copies an existing HTML dir
    # Creates new JS files in that dir's JS
    '''
    
    args = sys.argv
    if args[-1] != "1":
        help_str = "python3 make_BarSeqR_html_dir.py new_HTML_dir_path" + \
                   " inp_fit.tsv inp_t_score.tsv organism_name" + \
                   " exps_start_col 1\n" + \
                   "\nNote that exps_start_col refers to the column at" + \
                   " which the experiment values start, normally 5" 
        print(help_str)
    else:
        logging.basicConfig(level=logging.DEBUG)
        new_HTML_dir_path = args[1]
        inp_fit_score_fp = args[2]
        inp_t_score_fp = args[3]
        organism_name = args[4]
        col_start = args[5]
        create_new_HTML_dir(new_HTML_dir_path, 
                            inp_fit_score_fp, 
                            inp_t_score_fp, 
                            organism_name,
                            col_start)
        


    return None

if __name__ == "__main__":
    main()
