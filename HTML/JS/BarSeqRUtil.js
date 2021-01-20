


function BarSeqRComputeCorrelationBetweenChosenConditions() {
    // We use jstat to compute the correlation coefficient
    let num_selected = window.selected_conditions.length
    if (num_selected < 2) {
        alert("In order to compute correlation you need 2 selected conditions." 
                + " Currently " + num_selected.toString() + ". They are: " 
                + window.selected_conditions.join(", "))
    } else if (num_selected == 2) {

        // We get the matrix column index of the conditions
        let col_ix_a = window.exp2col['col2ix'][window.selected_conditions[0]]
        let col_ix_b = window.exp2col['col2ix'][window.selected_conditions[1]]
        
        // We create the two 'vectors' that hold the values
        let arr_a = []
        let arr_b = []
        for (let i=0; i<window.volcano_object['fit_mega_matrix'].length; i++) {
            arr_a.push(volcano_object['fit_mega_matrix'][i][col_ix_a])
            arr_b.push(volcano_object['fit_mega_matrix'][i][col_ix_b])
        }
        let corcoef = jStat.corrcoeff(arr_a, arr_b)
        return corcoef
    }
}

function BarSeqRComputeMeanof(typ) {
    // typ: str ('gene', 'cond' )
    // We use jStat lib to compute the mean
    
    if (typ == "gene") {
        let last_gene = window.last_plotted_gene
        let row_num = gene2row["gene_id4_2rownum"][last_gene][0]
        let values = volcano_object["fit_mega_matrix"][row_num]  
        mean = jStat.mean(values)
        return mean
    } else if (typ == "cond") {

        let num_selected = window.selected_conditions.length
        if (num_selected == 0) {
            alert("In order to compute the mean you need at least a single" 
                  + " selected condition.")
        } else {
            if (num_selected == 2) {
                alert("Computing the mean of the first of two selected conditions")
            }
            // We get the matrix column index of the conditions
            let col_ix = window.exp2col['col2ix'][window.selected_conditions[0]]

            // We create the vector that holds the values
            let arr_a = []
            for (let i=0; i<window.volcano_object['fit_mega_matrix'].length; i++) {
                arr_a.push(volcano_object['fit_mega_matrix'][i][col_ix])
            }
            let mean = jStat.mean(arr_a)
            return mean 
        }
    }
}


async function BarSeqR_PlotGeneAgainstAllConditions(gene_id4_str) {
    /* gene_id4_str: (str) A unique id for the gene with
     *  4 identifiers (does not include description)
    //
    //  Function process:
            gene_id4_str -> row in t-score matrix along with original
                        index of row
            row in t score matrix is sorted from lowest to highest 
            while original indices are maintained
            
            This list of lists 
            <<t score of condition, 1, index of condition>,
             <t score of condition, 2, index of condition>,
             ...,>
            is passed into the visualizer, where it is converted into 
            actual points on the plot.
    
    */

    // We note which gene 
    updateGeneInfoBox(gene_id4_str)
    window.last_plotted_gene = gene_id4_str

    let gene_v_cond_row_info = getRowInfoForGeneAgainstCondition(gene_id4_str)

    let gene_plot_object = {
        "min_x": 0.1,
        "max_x": gene_v_cond_row_info.length + 1,
        "min_y": gene_v_cond_row_info[0][1],
        "max_y": gene_v_cond_row_info[gene_v_cond_row_info.length - 1][1],
    }

    let new_ret_d = await refreshScatterPlotAxes('graph-svg',
                           gene_plot_object,
                           gene_plot_axes)



    populateSVGWithScatterPoints('graph-svg',
                                 gene_v_cond_row_info,
                                 window.ret_d['x_ticks_list'],
                                 window.ret_d['y_ticks_list'],
                                 window.ret_d['x_axis_start'],
                                 window.ret_d['y_axis_start'],
                                 window.ret_d['x_axis_len'],
                                 window.ret_d['y_axis_len'],
                                 point_contains_data = true, 
                                 point_click_function = BarSeqRShowInfoRelatedToPoint,
                                 point_radius = 7
                                )

    //Add Mean Fitness to SVG
    let mf_text = "Mean Fitness: " + (BarSeqRComputeMeanof('gene')).toString().slice(0,6)
    makeText(d3svg, 'normal', 20, 20, 20, mf_text, 'black', id_txt=null) 

    return null;

}

function updateGeneInfoBox(gene_id4_str) {
        
    /* gene_id4_str: (str) A unique id for the gene with
     *  4 identifiers (does not include description)
     */

    let gene_desc = gene2row["gene_id4_2rownum"][gene_id4_str][1];
    let gene_display_dobj = document.getElementById('selected-gene-display');
    gene_display_dobj.innerHTML = "Gene description: " + gene_desc;

    

}

function getRowInfoForGeneAgainstCondition(gene_id4_str) {

    /* gene_id4_str: (str) A unique id for the gene with
     *  4 identifiers (does not include description)
     */
    
    if (!(gene_id4_str in gene2row["gene_id4_2rownum"])) {
        throw gene_id4_str + " not found in gene ids object, cannot plot."
    }
    let row_num = gene2row["gene_id4_2rownum"][gene_id4_str][0]
    let fit_score_row = volcano_object["fit_mega_matrix"][row_num]  
    let fit_score_and_condition_ix_list = []
    for (let i =0; i<fit_score_row.length; i++) {
        fit_score_and_condition_ix_list.push([fit_score_row[i], i])
    }
    // Now we sort the above list by the t score value of the sublists
    fit_score_and_condition_ix_list.sort((a,b) => {
        if (a[0] >= b[0]) {
            return 1
        } else {
            return -1
        }
    })

    // Now we add the internal numbering for plotting's sake
    for (let i=0; i<fit_score_and_condition_ix_list.length; i++) {
        let point_data = {
            "typ": "gene_cond",
            "condition_column": fit_score_and_condition_ix_list[i][1],
            "y_val": fit_score_and_condition_ix_list[i][0]
        }
        fit_score_and_condition_ix_list[i] = [i,
                                              fit_score_and_condition_ix_list[i][0],
                                            point_data
                                            ]
    }
    return fit_score_and_condition_ix_list
}


// Here we do functions like populating the table 
// with the experiments

function BarSeqRPopulateTable(col2ix_d, ix2col_d) {
 /* 
  * We go from sorted alphabetical setName to
  *     number within list
  * col2ix_d:
  *     setName -> index
  * ix2col_d:
  *     index of column -> setName
  *
  */
   
    // We use the function from
    // LayoutUtil:
    // createTableWithRows(table_id, rows_info)
    let setnamelist = Object.keys(col2ix_d) 
    // we sort alphabetically
    setnamelist.sort()
    let rows_info = createRowsInfoFromIx2Col(setnamelist, col2ix_d); 
    //console.log(rows_info)
    createTableWithRows('exp-table', rows_info)

}


function createRowsInfoFromIx2Col(set_name_list, col2ix) {
     /*
    *
    * Args:
    *   set_name_list list<set_name (str)>
    *   col2ix: (object)
    *       column_name (str) -> index (str)
    *
    * Returns:
     *     rows_info: list<row_info>
     *          row_info: list<cell_info, cell_info, ...>
     *              cell_info: Object
     *                  type: (str) ["link" or "text"]
     *                  inner_text: (str) The text inside the cell
     *                  id: (str) The DOM id
     *                  [color]: Color of text
     *                  [textDecoration]: How the text will look
     *                  
     *                  IF type is link:
     *                      [href]: Link to go to
     *                      [func]: The let onclick function
     *                          [func_params]: If func, then func_params
     *                              must be present.
  */

    let rows_info = []

    for (let i = 0; i<set_name_list.length; i++) {

        let row_info = []
        /*
        let index_cell_info = {
            "type": "text",
            "tag_type": "p",
            "inner_text": index_list[i],
            "id": "col-num-" + index_list[i],
            "color": "blue"
        }
        */
        let set_name_info = {
            "type": "link",
            "inner_text": set_name_list[i] ,
            "id": "col-set-" + col2ix[set_name_list[i]],
            "color": "blue",
            "func": BarSeqRAddToSelected,
            "func_params":['chosen-display-div', set_name_list[i]] 
        }
        row_info.push(set_name_info)
        //row_info.push(index_cell_info)
        //console.log(row_info)
        rows_info.push(row_info)

    }

    return rows_info
}



function BarSeqRAddToSelected(inp_list) {
    /*
     * Args:
     *   inp_list: list<selected_div_id, selected_condition_str>
     *   selected_div_id: (str) Id for selected display div
     *   selected_condition_str: (str)
     */
    
    let selected_div_id = inp_list[0]
    let selected_condition_str = inp_list[1]
    //console.log(selected_div_id)
    //console.log(selected_condition_str)
   
    if (window.selected_conditions == null) {
        window.selected_conditions = []
    }
    //Note that the selected list is a global list
    let num_selected = window.selected_conditions.length
    if (num_selected >= 2) {
        alert("maximum number of selected conditions to compare is 2")
    } else {
        if (num_selected == 1) {
            if (selected_condition_str == window.selected_conditions[0]) {
                alert("Condition already selected")
                return;
            }
        }
        window.selected_conditions.push(selected_condition_str)

        let selected_div_dobj = document.getElementById(selected_div_id)
        while (selected_div_dobj.firstChild) {
            selected_div_dobj.removeChild(selected_div_dobj.firstChild);
        }

        //console.log(selected_div_dobj)
        let color_d = {
            "0": "#D2691E",
            "1": "#48D1CC"
        }
        for (let i = 0; i<window.selected_conditions.length; i++) {
            let current_p_tag = document.createElement("p")
            current_p_tag.innerHTML = window.selected_conditions[i]
            current_p_tag.style.color = color_d[i.toString()]
            selected_div_dobj.appendChild(current_p_tag)
        }

    }
}

function BarSeqRClearSelected() {

    // First we remove the actual selected
    window.selected_conditions = [] 
    // Then we remove the children in the display
    BarSeqRRemoveChildrenNodes('chosen-display-div')
    BarSeqRRemoveChildrenNodes('selected-point-display')
    BarSeqRRemoveChildrenNodes('selected-gene-display')

    //console.log("Cleared selected conditions.")
    refreshScatterPlotAxes('graph-svg',
                           volcano_object,
                           SVGGraphAxes)
}

function BarSeqRRemoveChildrenNodes(dobj_id) {

    // We remove the children 
    let selected_div_dobj = document.getElementById(dobj_id)
    while (selected_div_dobj.firstChild) {
        selected_div_dobj.removeChild(selected_div_dobj.firstChild);
    }
    return selected_div_dobj
}


function BarSeqRPrepareOnClicks() {
            let clr_btn = document.getElementById('clear-selected-btn-div')
            clr_btn.onclick = function() {
                BarSeqRClearSelected()
            }
            let volcano_btn = document.getElementById('volcano-btn-div')
            volcano_btn.onclick = function() {
                BarSeqRPlotFirstSelectedCondition()
            }
            let compare_plot_btn = document.getElementById('compare-btn-div')
            compare_plot_btn.onclick = function() {
                BarSeqRCompareConditions()
            }

    
    //deep copy of SVGGraphAxes for Compare Conditions Plot:

    if (!(window.hasOwnProperty('compare_axes_info'))) {
        window.compare_axes_info = JSON.parse(JSON.stringify(SVGGraphAxes));
        compare_axes_info["y_i"]["y_title_i"]["label"] = "Fitness Score"
        compare_axes_info["y_i"]["y_title_i"]["style_i"]["fontColor"] = "#48D1CC"
        compare_axes_info["x_i"]["x_title_i"]["style_i"]["fontColor"] = "#D2691E"
    }
    if (!(window.hasOwnProperty('gene_plot_axes'))) {
        window.gene_plot_axes = JSON.parse(JSON.stringify(SVGGraphAxes));
        gene_plot_axes["x_i"]["x_title_i"]["label"] = "Condition"
        gene_plot_axes["y_i"]["y_title_i"]["label"] = "Fitness Score"
        gene_plot_axes["y_i"]["y_title_i"]["style_i"]["fontColor"] = "black"
        gene_plot_axes["x_i"]["x_title_i"]["style_i"]["fontColor"] = "black"
    }
}

function BarSeqRgetConditionFitnessVsTScoreList(condition_str) {
    /*
     * We use the two matrices (fitness mega matrix)
     * and (t score mega matrix) to get a list of 
     * (<fitness (Num), t score (Num), point_data (Object)>)
     * for all columns for a condition.
     * 
     * point_data:
     *
    //      typ: (str) from ["volcano", "compare"] This will be 'volcano'
    //      row_num: (Num) int representing row from which
    //               point came
    //      x_val: Num
    //      y_val: Num
     *
     *
     * The matrices exist at window.volcano_d['t_mega_matrix']
     *      and window.volcano_d['fit_mega_matrix']
     *  The transition from condition string to index
     *  exists at window.exp2col['col2ix']
     *
     *
     *
     */ 

    let col_ix = window.exp2col['col2ix'][condition_str]

    let fit_to_t_list = []

    for (let i=0; i<window.volcano_object['fit_mega_matrix'].length; i++) {
        let point_data = {
            "typ": "volcano",
            "row_num": i.toString(),
            "x_val": window.volcano_object['fit_mega_matrix'][i][col_ix],
            "y_val": window.volcano_object['t_mega_matrix'][i][col_ix]
        }
        fit_to_t_list.push([window.volcano_object['fit_mega_matrix'][i][col_ix],
                            window.volcano_object['t_mega_matrix'][i][col_ix],
                            point_data])
    }

    return fit_to_t_list

}

async function BarSeqRPlotCondition(condition_str, clear_svg=true) {

    /*
    if (clear_svg) {
        await refreshScatterPlotAxes('graph-svg', 
                               volcano_object, 
                               SVGGraphAxes)
    }
    */

    // ret_d is attached to the windo
    let fit_to_t_list = await BarSeqRgetConditionFitnessVsTScoreList(condition_str)

    //console.log(window.ret_d)

    // function from 'ScatterPlot.js'
    populateSVGWithScatterPoints('graph-svg',
                                 fit_to_t_list,
                                 window.ret_d['x_ticks_list'],
                                 window.ret_d['y_ticks_list'],
                                 window.ret_d['x_axis_start'],
                                 window.ret_d['y_axis_start'],
                                 window.ret_d['x_axis_len'],
                                 window.ret_d['y_axis_len'],
                                 point_contains_data = true, 
                                 point_click_function = BarSeqRShowInfoRelatedToPoint
                                )

    //Add Mean Fitness to SVG
    let mf_text = "Mean Fitness: " + (BarSeqRComputeMeanof('cond')).toString().slice(0,6)
    makeText(d3svg, 'normal', 20, 20, 20, mf_text, 'black', id_txt=null) 

}

function BarSeqRPlotFirstSelectedCondition() {
    let condition_to_plot = null;
    if (window.selected_conditions.length == 0) {
        alert('No selected conditions to plot.')
        return
    } else if (window.selected_conditions.length == 2) {
        alert('Plotting the first condition selected out of 2')
    }
    refreshScatterPlotAxes('graph-svg',
                           volcano_object,
                           SVGGraphAxes)
    condition_to_plot = window.selected_conditions[0]
    BarSeqRPlotCondition(condition_to_plot)
}

function BarSeqRGetGeneInfoFromRow(row_num) {
    // The dict we'll use is window.gene2row['rownum2gene']

    let gene_info = window.gene2row['rownum2gene'][row_num.toString()]
    let gene_info_list = gene_info.split('|')

    return gene_info_list



}

function BarSeqRShowInfoRelatedToPoint(inp_d) {

    // inp_d: Object which contains the following keys: 
    //      typ: (str) from ["volcano", "compare", "gene_cond"] 
    //          IF typ volcano, compare:
    //              row_num: (Num) int representing row from which
    //                       point came
    //              x_val: Num
    //              y_val: Num
    //          IF typ gene_cond:
    //              [condition_column]: (Num) int representing
    //                      column from which point came
    //              [y_val]: Num (fitness value)

    let display_dobj = BarSeqRRemoveChildrenNodes('selected-point-display')
    let gene_info_list = null 
    if (inp_d['typ'] == "volcano" || inp_d['typ'] == "compare") {

        gene_info_list = BarSeqRGetGeneInfoFromRow(inp_d['row_num'])

        gene_display = document.createElement("p")
        gene_plot_link = document.createElement("a")

        gene_display.innerHTML = "Gene Description: " + gene_info_list[4] +
                                ". Locus Tag: " + gene_info_list[1] + 
                                ". Gene SysName: " + gene_info_list[2] + "."
        gene_plot_link.innerHTML = "Plot Gene"
        gene_plot_link.style.textDecoration = "underline"
        gene_plot_link.style.color = "blue"
        gene_plot_link.style.cursor = "pointer"
        gene_plot_link.onclick = function () {
            gene_id4_str = gene_info_list.slice(0,4).join('|')
            console.log(gene_id4_str)
            BarSeqR_PlotGeneAgainstAllConditions(gene_id4_str)
        }
        display_dobj.appendChild(gene_display)
        display_dobj.appendChild(gene_plot_link)

        if (inp_d['typ'] == "volcano") {
            fit_display = document.createElement("p")
            t_display = document.createElement("p")
            fit_display.innerHTML = "Fitness: " + inp_d['x_val'].toString()
            t_display.innerHTML = "T score: " + inp_d['y_val'].toString()
            display_dobj.appendChild(fit_display)
            display_dobj.appendChild(t_display)
        } else {
            x_fit_display = document.createElement("p")
            y_fit_display = document.createElement("p")
            x_t_display = document.createElement("p")
            y_t_display = document.createElement("p")
            x_fit_display.innerHTML = "X Fitness: " + inp_d['x_val'].toString()
            x_t_display.innerHTML = "X T-Score (abs): " + inp_d['x_t_score'].toString()
            y_fit_display.innerHTML = "Y Fitness: " + inp_d['y_val'].toString()
            y_t_display.innerHTML = "Y T-Score (abs): " + inp_d['y_t_score'].toString()
            display_dobj.appendChild(x_fit_display)
            display_dobj.appendChild(x_t_display)
            display_dobj.appendChild(y_fit_display)
            display_dobj.appendChild(y_t_display)
        }

    } else if (inp_d['typ'] == "gene_cond") {
        let condition_link = document.createElement("a")
        let col_str = exp2col['ix2col'][inp_d['condition_column']]
        condition_link.innerHTML = col_str 
        condition_link.onclick = function() {
            BarSeqRAddToSelected(['chosen-display-div', col_str])
        }
        let y_fit_display = document.createElement("p")
        y_fit_display.innerHTML = "Fitness Value: " + inp_d['y_val'].toString()
        display_dobj.appendChild(condition_link)
        display_dobj.appendChild(y_fit_display)
    } else {
        console.log("Cannot recognize point type - no display")
    }

}

function BarSeqRCompareConditions() {


    let sC = window.selected_conditions
    if ( sC.length < 2) {
        alert("You must select two conditions to compare")
        return;
    } else {
        refreshScatterPlotAxes('graph-svg',
                           volcano_object,
                           SVGGraphAxes)
        BarSeqRPlotComparedConditions(sC[0], sC[1])

        //Add Correlation Fitness to SVG
        let fc_text = "Fitness Correlation: " + 
            BarSeqRComputeCorrelationBetweenChosenConditions().toString().slice(0,6);
        makeText(d3svg, 'normal', 20, 20, 20, fc_text, 'black', id_txt=null) 
    }
}

async function BarSeqRPlotComparedConditions(cond_a, cond_b) {
    /*
     * Args:
     *      cond_a, cond_b: (str) Condition values
     *
     */

    // We get the matrix column index of the conditions
    let col_ix_a = window.exp2col['col2ix'][cond_a]
    let col_ix_b = window.exp2col['col2ix'][cond_b]

    let condition_info_obj = BarSeqRGetConditionListAndMinMax(col_ix_a, col_ix_b)

    //console.log(condition_info_obj)

    
    let new_ret_d = await refreshScatterPlotAxes('graph-svg',
                           condition_info_obj,
                           compare_axes_info)
    
    
    
    // function from 'ScatterPlot.js'
    populateSVGWithScatterPoints('graph-svg',
                                 condition_info_obj['point_list'],
                                 new_ret_d['x_ticks_list'],
                                 new_ret_d['y_ticks_list'],
                                 new_ret_d['x_axis_start'],
                                 new_ret_d['y_axis_start'],
                                 new_ret_d['x_axis_len'],
                                 new_ret_d['y_axis_len'],
                                 point_contains_data = true, 
                                 point_click_function = BarSeqRShowInfoRelatedToPoint
                                )
    

    
}

function BarSeqRGetConditionListAndMinMax(col_ix_a, col_ix_b) {
    // We get the list of values to plot and min and max
    // Args: indeces within fitness matrix of columns for conditions
    //
    // Returns:
    //     Object like volcano_object. Keys listed at bototm of func
    //
    

    let min_x = 0
    let min_y = 0
    let max_x = 0
    let max_y = 0

    let fita2fitb_list = [];

    for (let i=0; i<window.volcano_object['fit_mega_matrix'].length; i++) {

        let x_val = window.volcano_object['fit_mega_matrix'][i][col_ix_a]
        let y_val = window.volcano_object['fit_mega_matrix'][i][col_ix_b]
        let x_t_score = window.volcano_object['t_mega_matrix'][i][col_ix_a]
        let y_t_score = window.volcano_object['t_mega_matrix'][i][col_ix_b]

        // Checking for x and y min/max
        if (x_val < min_x) {
            min_x = x_val
        } else if (x_val > max_x) {
            max_x = x_val
        }

        if (y_val < min_y) {
            min_y = y_val
        } else if (y_val > max_y) {
            max_y = y_val
        }

        let point_data = {
            "typ": "compare",
            "row_num": i.toString(),
            "x_val": x_val,
            "y_val": y_val,
            "x_t_score": x_t_score,
            "y_t_score": y_t_score
        }
        fita2fitb_list.push([x_val,
                            y_val,
                            point_data])
    }

    /*
    min_x = Math.ceil(min_x - 1)
    max_x = Math.floor(max_x + 1)
    min_y = Math.ceil(min_y - 1)
    max_y = Math.floor(max_y + 1)
    */

    return {
        "point_list": fita2fitb_list,
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y
    }


}


