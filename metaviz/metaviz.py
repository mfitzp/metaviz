# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
from copy import copy
from optparse import Values, OptionParser
from collections import defaultdict
from biocyc import biocyc, Pathway, Protein, Reaction, Compound, Gene
import pydot

import logging

PRUNE_ALL = lambda a, b, c, d: (a, b, c)
PRUNE_IDENTICAL = lambda a, b, c, d: (a, b, c, d)

# Internal URLS
COMPOUND_URL = 'pathomx://db/compound/%s/view'
PATHWAY_URL = 'pathomx://db/pathway/%s/view'
REACTION_URL = 'pathomx://db/reaction/%s/view'
PROTEIN_URL = 'pathomx://db/protein/%s/view'
GENE_URL = 'pathomx://db/gene/%s/view'

# Paper sizes for print scaling printing
METAPATH_PAPER_SIZES = {
    'None': (-1, -1),
    'A0': (33.11, 46.81),
    'A1': (23.39, 33.11),
    'A2': (16.54, 23.39),
    'A3': (11.69, 16.54),
    'A4': (8.27, 11.69),
    'A5': (5.83, 8.27),
}

REACTION_DIRECTIONS = {
    'LEFT-TO-RIGHT': 'forward',
    'RIGHT-TO-LEFT': 'back',
    'REVERSIBLE': 'both',
    'IRREVERSIBLE-LEFT-TO-RIGHT': 'forward',
    'IRREVERSIBLE-RIGHT-TO-LEFT': 'back',
    'PHYSIOL-LEFT-TO-RIGHT': 'forward',
    'PHYSIOL-RIGHT-TO-LEFT': 'back'
    }

SECONDARY_METABOLITES = [
    # Nuceleosides
    'AMP', 'ADP', 'ATP',
    'CMP', 'CDP', 'CTP',
    'GMP', 'GDP', 'GTP',
    'UMP', 'UDP', 'UTP',
    'TMP', 'TDP', 'TTP',
    # Deoxy-nucleosides (only
    'DAMP', 'DADP', 'DATP',
    'DCMP', 'DCDP', 'DCTP',
    'DGMP', 'DGDP', 'DGTP',
    'DUMP', 'DUDP', 'DUTP',
    'DTMP', 'DTDP', 'DTTP',
    # Reducing
    'NAD', 'NAD-P-OR-NOP', 'NADP',
    'NADH', 'NADH-P-OR-NOP', 'NADPH',
    'CPD-653', 'CPD0-2472',  # NADHX-S and -R
    'FAD',
    'FADH', 'FADH2',
    # Protons
    'PROTON', 'Donor-H2', 'Acceptor',
    # Molecules
    'CARBON-DIOXIDE', 'WATER', 'OXYGEN-MOLECULE',
    # Metal ions
    'CA', 'CA+2',
    'FE', 'FE+2', 'FE+3',
    #Inorganic phosphate
    'Pi', 'PPI', 'P3I', 'P4I',
    # Miscellaneous
    'HCO3',  # SO3, #AMMONIA
    'Menaquinones',
    'Unspecified-Degradation-Products', 'Demethylated-methyl-acceptors',
    # Ubiquino
    'Ubiquinones', 'Ubiquinols',  # 'UBIQUINONE-8',#'CPD-9956','CPD-9958',
    # Co-Enzymes
    'CO-A', 'BIOTIN',
    #'BIOTIN','THIAMINE-PYROPHOSPHATE', 'PYRIDOXAL_PHOSPHATE', 'THF',

    # Energy
    'UV-Light', 'Light',
    ]


class ReactionIntermediate(object):
    pass

class DummyReaction(object):
    pass

class DummyPathway(object):
    def __str__(self):
        if self.name:
            return self.name
        else:
            return self.id

def add_clusternodes(clusternodes, cluster_key, keys, nodes):
    for key in keys:
        clusternodes[cluster_key][key].extend(nodes)
    return clusternodes


def get_compound_color(analysis, m):
    return analysis[m.id][0] if m.id in analysis else '#cccccc'


def get_pathway_color(analysis, m):
    return analysis[m.id][0] if m.id in analysis else '#cccccc'


def get_reaction_color(analysis, r):
    # For reactions, we need gene and protein data (where it exists)
    colors = []
    if hasattr(r, 'enzymes'):  # Dummy ReactionIntermediates will not; til that's fixed!
        for p in r.enzymes:
            if p in analysis:
                colors.append(analysis[p][0])

            if hasattr(p, 'genes'):
                for g in p.genes:
                    if g in analysis:
                        colors.append(analysis[g][0])

    if colors == []:
        colors = ['#cccccc']  # Mid-grey      

    return '"%s"' % ':'.join(colors)

def major_compounds(l):
    major = [o for o in l if o.id not in SECONDARY_METABOLITES]
    return major

def minor_compounds(l):
    minor = [o for o in l if o.id in SECONDARY_METABOLITES]
    return minor

def get_filename_with_counter(filename):
    fn, ext = os.path.splitext(filename)
    return fn + "-%s" + ext

def remove_html_markup(s):
    tag = False
    quote = False
    out = ""

    for c in s:
        if c == '<' and not quote:
            tag = True
        elif c == '>' and not quote:
            tag = False
        elif (c == '"' or c == "'") and tag:
            quote = not quote
        elif not tag:
            out = out + c

    return out

def generate(pathways, reactions=[], analysis=None, organism='HUMAN',
        cluster_by='pathway',
        
        show_enzymes=True,
        show_secondary=True,
        show_molecular=True,
        show_network_analysis=True,
        show_gibbs=True,

        highlightpathways=True,
        highlightregions=True,

        show_pathway_links=False):
    
    # Initialise the BioCyc database accessor
    biocyc.set_organism(organism)

    # Convert pathway list to list of BioCyc pathway objects (unless they already are)
    pathways = [ p if type(p) == Pathway else biocyc.get( str(p) ) for p in pathways ]
    
    if analysis:
        # Convert analysis data to BioCyc pathway objects
        updated_analysis = {}
        for k,v in analysis.items():
            if type(k) in [Gene, Compound, Protein]:
                updated_analysis[k] = v
            else:
                updated_analysis[ biocyc.get( str(k) ) ] = v
            
        analysis = updated_analysis

    # Wrap solo reactions in fake meta-pathways
    reactions = [ r if type(r) == Reaction else biocyc.get( str(r) ) for r in reactions ]
    for r in reactions:
        p = DummyPathway()
        p.id = r.id
        p.name = r.name
        p.reactions = [r]
        p.compounds = r.compounds
        pathways.append(p)
        

    # Internode counter (create dummy nodes for split compounds)
    intno = 0
    visible = True
    splines = True

    # Pathway colour list
    #                   black       blue        green       red         orange      purple      yellow      pink
    colors = ['#000000', '#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#392c85', '#ffff00', '#ff4398']  # sat 100
    colorslight = ['#aaaaaa', '#ccccff', '#b3ffb3', '#ffb3b3', '#ffedcc', '#ffb3fe', '#ffffcc', '#ffcce3']  # sat 20
    colorslighter = ['#cccccc', '#f2f2ff', '#f2fff2', '#fff2f2', '#fffbf2', '#fff2ff', '#fffff2', '#fff2f8']  # sat 5

    prunekey = PRUNE_IDENTICAL if (show_enzymes or show_secondary) else PRUNE_ALL

    # Subgraphs of metabolic pathways
    nodes = list()
    edges = list()
    edgesprune = list()

    inter_node = 0
    itr = 0

    clusternodes = dict()
    clusternodes['pathway'] = defaultdict(list)
    clusternodes['compartment'] = defaultdict(list)

    edgecluster = dict()
    edgecluster['pathway'] = defaultdict(list)
    edgecluster['compartment'] = defaultdict(list)

    clusters = dict()
    clusters['pathway'] = set(pathways)
    clusters['compartment'] = set()

    cluster_key = cluster_by
    clusters['compartment'].add('Non-compartmental')  # Need to override the color on this later

    # Store alternative pathways for reactions, for use when pruning deletes reactions
    pathway_edges_alternates = defaultdict(tuple)  # Dict of tuples pathway => 

    for p in pathways:
        for r in p.reactions:
            if r is None:
                continue
        # Check that this edge is between items in one of the specified pathways
            # FIXME: biocyc compartment support
            #compartments = [c for pr in r.enzymes for c in pr.compartments]
            compartments = []
            if compartments == []:
                compartments = ['Non-compartmental']
            clusters['compartment'] |= set(compartments)  # Add to cluster set
            # Store edge cluster data (reaction)
            edgecluster['pathway'][r].append(p)
            edgecluster['compartment'][r].extend(compartments)

            nmtins = set()
            nmtouts = set()
            
            
            if not hasattr(r,'direction'):
                r.direction = 'REVERSIBLE'

            for mtin in major_compounds(r.compounds_left):
                for mtout in major_compounds(r.compounds_right):
                    # Use a in/out/enzyme tuple to delete duplicates
                    if prunekey(mtin, mtout, r.direction, r.enzymes) in edgesprune:
                        continue
                    else:
                        edgesprune.append(prunekey(mtin, mtout, r.direction, r.enzymes))
                    nmtins.add(mtin)
                    nmtouts.add(mtout)

            if nmtins and nmtouts:

                mtins = list(nmtins)
                mtouts = list(nmtouts)

                # Make a copy of the reaction object, so we can add link data
                inter_react = DummyReaction()
                inter_react.id = r.id
                inter_react.name = ''
                inter_react.enzymes = r.enzymes
                inter_react.type = 'dummy'
                inter_react.direction = r.direction

                edgecluster['pathway'][inter_react].append(p)
                edgecluster['compartment'][inter_react].extend(compartments)

                # If multiple ins/outs create dummy split-nodes
                # RXNINXX, RXNOUTXX
                if len(mtins) > 1:
                    intno += 1  # Increment no
                    inter_node = ReactionIntermediate()
                    inter_node.id = "DUMMYRXN-IN%d" % intno
                    inter_node.type = 'dummy'
                    
                    for mtin in mtins:
                        edges.append([inter_react, mtin, inter_node, visible])
                        clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [mtin])
                        clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [mtin])

                    clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [inter_node])
                    clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [inter_node])

                    nodes.append([inter_node, False, visible])
                    # Overwrite with the dummy name, use this as the basis of the main detail below
                    mtin = inter_node
                else:
                    mtin = mtins[0]

                if len(mtouts) > 1:
                    intno += 1  # Increment no
                    inter_node = ReactionIntermediate()
                    inter_node.id = "DUMMYRXN-OUT%d" % intno
                    inter_node.type = 'dummy'
                    
                    for mtout in mtouts:
                        edges.append([inter_react, inter_node, mtout, visible])
                        clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [mtout])
                        clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [mtout])

                    clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [inter_node])
                    clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [inter_node])

                    nodes.append([inter_node, False, visible])
                    # Overwrite with the dummy name, use this as the basis of the main detail below
                    mtout = inter_node
                else:
                    mtout = mtouts[0]

                edges.append([r, mtin, mtout, visible])

                # Store clustering data for layout
                clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [mtin])
                clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [mtin])
                clusternodes = add_clusternodes(clusternodes, 'pathway', [p], [mtout])
                clusternodes = add_clusternodes(clusternodes, 'compartment', compartments, [mtout])

    pathway_compounds = [major_compounds(p.compounds) for p in pathways] # Get major only for the nodes
    pathway_compounds = [item for sublist in pathway_compounds for item in sublist]

    # id,type,names
    for m in pathway_compounds:

        # It's in one of our pathways (union)
        fillcolor = False

        if analysis and m in analysis:
            # We found it by one of the names
            fillcolor = analysis[m]

        # This node is in one of our pathways, store it
        nodes.append([m, fillcolor, visible])

    # Add pathway annotations
    if show_pathway_links:

        visible_reactions = [r for r, x1, x2, x3 in edges]
        visible_nodes = [n for n, x1, x2 in nodes]

        pathway_annotate = set()
        pathway_annotate_dupcheck = set()
        
        pathway_reactions = [p.reactions for p in pathways]
        pathway_reactions = [item for sublist in pathway_reactions for item in sublist]
        
        for r in pathway_reactions:

        # Check that a reaction for this isn't already on the map
            if r not in visible_reactions:
                # Now find out which end of it is (one side's compounds [or both])
                for p in r.pathways:
                    pathway_node = ReactionIntermediate(**{'id': '%s' % p.id, 'name': p.name, 'type': 'pathway'})

                    for mt in r.mtins:
                        if mt in visible_nodes and (p, mt) not in pathway_annotate_dupcheck:  # Compound is already on the graph
                            mp = mt.pathways[0]
                            pathway_annotate.add((p, mp, pathway_node, mt, pathway_node, r.dir))
                            pathway_annotate_dupcheck.add((p, mt))
                            break

                    for mt in r.mtouts:
                        if mt in visible_nodes and (p, mt) not in pathway_annotate_dupcheck:  # Compound is already on the graph
                            mp = mt.pathways[0]
                            pathway_annotate.add((p, mp, pathway_node, pathway_node, mt, r.dir))
                            pathway_annotate_dupcheck.add((p, mt))
                            break

    
        for p, mp, pathway_node, mtin, mtout, dir in list(pathway_annotate):
            itr += 1
            #nodepathway[mp].append(pathway_node)
            inter_react = ReactionIntermediate(**{'id': "DUMMYPATHWAYLINK-%s" % itr, 'type': 'dummy', 'dir': dir, 'pathways': [mp]})
            edges.append([inter_react, mtin, mtout, True])

            if analysis:  # and mining:
                # Not actually used for color, this is a ranking value (bud-sized on pathway link)
                p_compound_scores = [analysis[c] for c in p.compounds if c in analysis]
                if p_compound_scores:
                    fillcolor = sum(p_compound_scores) / len(p_compound_scores)
                else:
                    fillcolor = None
                # fillcolor = max(1, 11-analysis['mining_ranked_remaining_pathways'].index( p.id ) ) if p.id in analysis['mining_ranked_remaining_pathways'] else 1
            else:
                fillcolor = None

            nodes.append([pathway_node, fillcolor, True])

    # Generate the analysis graph from datasets
    graph = pydot.Dot('\u200C', graph_type='digraph', sep="+15,+10", esep="+5,+5", labelfloat='false', outputMode='edgeslast', fontname='Calibri', splines=splines, gcolor='white', pad=0.5, model='mds', overlap="vpsc")  # , model='mds') #, overlap='ipsep', mode='ipsep', model='mds')
    subgraphs = list()
    clusterclu = dict()

    nodes_added = set()  # Store nodes that are added, can use simplified adding for subsequent pathways

    # Arrange layout grouping (e.g. by pathway, compartment, etc.)
    for sgno, cluster in enumerate(clusters[cluster_key]):
        clusterclu[cluster] = (sgno % 7) + 1

        if highlightregions:
            if cluster_by == 'compartment':
                bcolor = colorslight[clusterclu[cluster]]
                bgcolor = colorslighter[clusterclu[cluster]]
                style = 'rounded'
            else:
                bcolor = '#eeeeee'
                bgcolor = 'transparent'
                style = 'solid'
        else:
            bcolor = 'transparent'
            bgcolor = 'transparent'
            style = 'solid'

        subgraph = pydot.Cluster(str(sgno), label='%s' % cluster, graph_type='digraph', fontname='Calibri', splines=splines, color=bcolor, bgcolor=bgcolor, style=style, fontcolor=bcolor, labeljust='left', pad=0.5, margin=12, labeltooltip='%s' % cluster, URL='non')  # PATHWAY_URL % cluster.id )
        # Read node file of compounds to show
        # TODO: Filter this by the option specification
        for n in clusternodes[cluster_key][cluster]:
            subgraph.add_node(pydot.Node(n.id))
        graph.add_subgraph(subgraph)

    # Add nodes to map
    for m, node_color, visible in nodes:
        if m in nodes_added:  # Previously added, another pathway: use simplified add (speed up)
            graph.add_node(pydot.Node(m.id))
            continue  # Next

        label = ' '
        color = 'black'
        shape = 'rect'
        fontcolor = 'black'
        fillcolor = '#eeeeee'
        colorscheme = 'rdbu9'
        url = COMPOUND_URL
        width, height = 0.75, 0.5

                    
        if visible:
            style = 'filled'
        else:
            style = 'invis'

        if m.type == 'dummy':
            shape = 'point'
            fillcolor = 'black'
            border = 0
            width, height = 0.01, 0.01
            url = 'pathomx://null/%s'  # Null, don't navigate FIXME

        elif m.type == 'pathway':
            shape = 'point'
            label = '%s' % remove_html_markup(m.name)
            size = len(m._compounds)
            width, height = size / 24., size / 24.
            if node_color is False:
                fillcolor = '#cccccc'
                color = '#cccccc'
            else:
                color = node_color[0]
            border = 0
            url = PATHWAY_URL

        else:
            label = label = "%s" % remove_html_markup(m.name)  # {%s |{ | | } } "
            if node_color == False:
                if analysis:  # Showing data
                    fillcolor = '#ffffff'
                else:
                    fillcolor = '#cccccc'
            else:
                shape = 'box'
                style = 'filled'
                fontcolor = node_color[1]
                fillcolor = node_color[0]

            if show_network_analysis:
                border = min(len(m._reactions) / 2, 5)
            else:
                border = 0

        if show_molecular and hasattr(m, 'image'):
            label = ' '
            if analysis and node_color and isinstance(node_color[2], int):
                image = m.imagecolor % int(node_color[2])
            else:
                image = m.image
            style = 'solid'
            shape = 'none'
        else:
            image = False

        if image:
            graph.add_node(pydot.Node(m.id, width=width, height=height, image=image, style=style, shape=shape, color=color, penwidth=border, fontname='Calibri', colorscheme=colorscheme, fontcolor=fontcolor, fillcolor=fillcolor, label=label, labeltooltip=label, URL=url % m.id))  # http://metacyc.org/META/substring-search?object=%s
        else:
            graph.add_node(pydot.Node(m.id, width=width, height=height, style=style, shape=shape, color=color, penwidth=border, fontname='Calibri', colorscheme=colorscheme, fontcolor=fontcolor, fillcolor=fillcolor, label=label, labeltooltip=label, URL=url % m.id))  # http://metacyc.org/META/substring-search?object=%s

        nodes_added.add(m)
        
        
    # Add graph edges to the map
    style = ' '
    for r, origin, dest, visible in edges:
        label = list()
        arrowhead = 'normal'
        arrowtail = 'empty'
        color = '#888888'
        url = REACTION_URL
        length = 2
        penwidth = 1
        weight = 1
        dir = REACTION_DIRECTIONS[r.direction]

        # End of any edge touching a DUMMY-RXN is left blank
        if dest.type == 'dummy':
            arrowhead = 'none'
            length = 1.3

        if origin.type == 'dummy':
            arrowtail = 'none'
            length = 1.3

        if visible:
            style = ' '
        else:
            style = 'invis'

        if analysis:
            color = get_reaction_color(analysis, r)
            colorscheme = 'rdbu9'

        elif highlightpathways:
            #color=1+( ] % 11) # Length of colorscheme -1
            r_clusterclu = list(set(edgecluster[cluster_key][r]) & set(clusterclu))
            color = '"%s"' % ':'.join([colors[n] for n in sorted([clusterclu[c] for c in r_clusterclu])])
            colorscheme = 'paired12'

        if r.type != 'dummy' and show_enzymes:
            label.append('%s' % r.name)

            if hasattr(r, 'enzymes') and r.enzymes:
                if analysis:
                    prgenestr = ''
                    if hasattr(r,'enzymes'):
                        for pr in r.enzymes:
                            if pr in analysis:
                                prgenestr += '<font color="%s">&#x25C6;</font>' % analysis[pr][0]
                                if hasattr(pr, 'genes'):
                                    for g in pr.genes:
                                        if g in analysis:
                                            prgenestr += '<font color="%s">&#x25cf;</font>' % analysis[g][0]
                    label.append('%s' % prgenestr)  # pr.genes

            mcl = minor_compounds( r.compounds_left )
            mcr = minor_compounds( r.compounds_right )
            if show_secondary and len( mcl + mcr ) > 0:
                # Process to add colors if compound in db
                smtins, smtouts = [], []
                for sm in mcl:
                    if analysis and sm in analysis:
                        smtins.append('<font color="%s">%s</font>' % (analysis[sm][0], sm))  # We found it by one of the names
                    else:
                        smtins.append('%s' % sm)

                for sm in mcr:
                    if analysis and sm in analysis:
                        smtouts.append('<font color="%s">%s</font>' % (analysis[sm][0], sm))  # We found it by one of the names
                    else:
                        smtouts.append('%s' % sm)
                
                label.append('%s → %s' % (', '.join(smtins), ', '.join(smtouts)))

        #if show_network_analysis:
        #    width = min( len( r._pathways ), 5)
        #else:
        #    width = 1

        if show_gibbs and hasattr(r, 'gibbs'):
            penwidth = abs(r.gibbs['deltaG_w'])

        e = pydot.Edge(origin.id, dest.id, weight=weight, len=length, penwidth=penwidth, dir=dir, label='<' + '<br />'.join(label) + '>', colorscheme=colorscheme, color=color, fontcolor='#888888', fontsize='10', arrowhead=arrowhead, arrowtail=arrowtail, style=style, fontname='Calibri', URL=url % r.id, labeltooltip=' ')
        graph.add_edge(e)

    return graph


