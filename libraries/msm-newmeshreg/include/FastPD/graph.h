//#############################################################################
//#
//# FastPD Optimization (c) 2008.
//#
//# Max-flow Computation.
//#
//# Original Author: Vladimir Kolmogorov
//# Modification: Nikos Komodakis
//#
//# This is an implementation of the modified max-flow computation described in:
//# N. Komodakis, G. Tziritas, N. Paragios,
//# Fast, Approximately Optimal Solutions for Single and Dynamic MRFs,
//# IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2007.
//#
//# THE WORK IS ONLY FOR RESEARCH AND NON-COMMERCIAL PURPOSES. THE OPTIMIZATION
//# CODE IS PROTECTED FROM SEVERAL US/EU/CHINA PENDING PATENT APPLICATIONS.
//# IF YOU WOULD LIKE TO USE THIS SOFTWARE FOR COMMERCIAL PURPOSES OR LICENSING
//# THE TECHNOLOGY, PLEASE CONTACT:
//# ECOLE CENTRALE DE PARIS, PROF. NIKOS PARAGIOS (nikos.paragios@ecp.fr).
//#
//# If you intend to use this code or results obtained with it, the above
//# mentioned paper should be cited within your publication.
//#
//#############################################################################

//#############################################################################
//# Header for the original max-flow computation by Yuri Boykov and
//# Vladimir Kolmogorov.
//#############################################################################
/* graph.h */
/* Vladimir Kolmogorov (vnk@cs.cornell.edu), 2001. */

/*
	This software library is a modification of the maxflow algorithm
	described in

	An Experimental Comparison of Min-Cut/Max-Flow Algorithms
	for Energy Minimization in Computer Vision.
	Yuri Boykov and Vladimir Kolmogorov.
	In Third International Workshop on Energy Minimization
	Methods in Computer Vision and Pattern Recognition, September 2001

	This algorithm was originally developed at Siemens.
	The main modification is that two trees are used for finding
	augmenting paths - one grows from the source and the other
	from the sink. (The original algorithm used only the former one).
	Details will be described in my PhD thesis.

	This implementation uses an adjacency list graph representation.
	Memory allocation:
		Nodes: 22 bytes + one field to hold a residual capacity
		       of t-links (by default it is 'short' - 2 bytes)
		Arcs: 12 bytes + one field to hold a residual capacity
		      (by default it is 'short' - 2 bytes)
	(Note that arcs are always added in pairs - in forward and reverse directions)

	Example usage (computes a maxflow on the following graph):

		        SOURCE
		       /       \
		     1/         \2
		     /      3    \
		   node0 -----> node1
		     |   <-----   |
		     |      4     |
		     \            /
		     5\          /6
		       \        /
		          SINK

	///////////////////////////////////////////////////

	#include <stdio.h>
	#include "graph.h"

	void main()
	{
		Graph::node_id nodes[2];
		Graph *g = new Graph();

		nodes[0] = g -> add_node();
		nodes[1] = g -> add_node();
		g -> set_tweights(nodes[0], 1, 5);
		g -> set_tweights(nodes[1], 2, 6);
		g -> add_edge(nodes[0], nodes[1], 3, 4);

		Graph::flowtype flow = g -> maxflow();

		printf("Flow = %d\n", flow);
		printf("Minimum cut:\n");
		if (g->what_segment(nodes[0]) == Graph::SOURCE)
			printf("node0 is in the SOURCE set\n");
		else
			printf("node0 is in the SINK set\n");
		if (g->what_segment(nodes[1]) == Graph::SOURCE)
			printf("node1 is in the SOURCE set\n");
		else
			printf("node1 is in the SINK set\n");

		delete g;
	}

	///////////////////////////////////////////////////
*/

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "block.h"

/*
	Nodes, arcs and pointers to nodes are
	added in blocks for memory and time efficiency.
	Below are numbers of items in blocks
*/
#define NODE_BLOCK_SIZE 512
#define ARC_BLOCK_SIZE 1024
#define NODEPTR_BLOCK_SIZE 128

/*
	special constants for node->parent
*/
#define TERMINAL ( (Graph::arc *) 1 )    /* to terminal */
#define ORPHAN   ( (Graph::arc *) 2 )    /* orphan */

#define INFINITE_D 1000000000    /* infinite distance to the terminal */

class Graph {
public:

#define _MANY_LABELS_
#ifdef _MANY_LABELS_
    typedef int Label;
#else
    typedef unsigned char Label;
#endif

    typedef enum {
        SOURCE = 0,
        SINK = 1
    } termtype; /* terminals */

/* Type of edge weights.
   Can be changed to char, int, float, double, ... */
    typedef double captype;
    typedef captype Real;

/* Type of total flow */
    typedef float flowtype;

    typedef void *node_id;

/* interface functions */

/* Destructor */
    ~Graph() {}

/* Adds a node to the graph */
    void add_nodes() {
        for (int i = 0; i < _num_nodes; i++) {
            _nodes[i].first = NULL;
            _nodes[i].tr_cap = 0;
#ifndef _METRIC_DISTANCE_
            _nodes[i].conflict_time = -1;
#endif
        }
    }

/* Adds a bidirectional edge between 'from' and 'to'
   with the weights 'cap' and 'rev_cap' */
    void add_edges(int *pairs, int numpairs) {
        captype cap = 0;
        captype rev_cap = 0;

        // numpairs iterations (note i+=2)
        for (int i = 0; i < 2 * numpairs; i += 2) {
            // retrieval of the node_id's (node*)
            node_id from = &_nodes[pairs[i]];
            node_id to = &_nodes[pairs[i + 1]];

            // there are always two arcs inserted into the graph
            // 'from -> to' as well as 'to -> from'
            arc *a, *a_rev;

            // retrieve the two arcs from the _arcs array
            // the array is indexed in the same way as pairs
            // arcs[i]   = node[pairs[i]] -> node[pairs[i+1]]
            // arcs[i+1] = node[pairs[i+1]] -> node[pairs[i]]
            a = &_arcs[i];
            a_rev = a + 1;

            // sister holds the reverse arcs
            a->sister = a_rev;
            a_rev->sister = a;

            // next points to the next outgoing arc/edge in from
            // since there is exactly one outgoing edge, it points
            // to itself!
            a->next = ((node *) from)->first; // i.e. a->next = a
            ((node *) from)->first = a;

            // here holds the same as above
            a_rev->next = ((node *) to)->first;
            ((node *) to)->first = a_rev;

            // head is pointing to the destination node
            a->head = (node *) to;
            a_rev->head = (node *) from;

            // and finally, we configure the capacities
            a->r_cap = cap;
            a_rev->r_cap = rev_cap;
        }
    }

/* Sets the weights of the edges 'SOURCE->i' and 'i->SINK'
   Can be called at most once for each node before any call to 'add_tweights'.
   Weights can be negative */
    void set_tweights(node_id i, captype cap_source, captype cap_sink) {
        flow += (cap_source < cap_sink) ? cap_source : cap_sink;
        ((node *) i)->tr_cap = cap_source - cap_sink;
    }

/* Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights
   Can be called multiple times for each node.
   Weights can be negative */
    void add_tweights(node_id i, captype cap_source, captype cap_sink) {
        captype delta = ((node *) i)->tr_cap;
        if (delta > 0) cap_source += delta;
        else cap_sink -= delta;
        flow += (cap_source < cap_sink) ? cap_source : cap_sink;
        ((node *) i)->tr_cap = cap_source - cap_sink;
    }


    /* After the maxflow is computed, this function returns to which
	   segment the node 'i' belongs (Graph::SOURCE or Graph::SINK) */
    termtype what_segment(node_id i);

/* Computes the maxflow. Can be called only once. */
    flowtype apply_maxflow(int init_on) {
        node *i, *j, *current_node = NULL;
        arc *a;
        nodeptr *np, *np_next;

        if (init_on)
            maxflow_init();
        nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);

//int num_iter = 0;					//ADDED BY BEN
//while ( num_iter < 10000 )		//ADDED BY BEN

        while (1) {
            //num_iter++;					//ADDED BY BEN

            if (i = current_node) {
                i->next = NULL; /* remove active flag */
                if (!i->parent) i = NULL;
            }
            if (!i) {
                if (!(i = next_active())) break;
            }

            /* growth */
            if (!i->is_sink) {
                /* grow source tree */
                for (a = i->first; a; a = a->next)
                    if (a->r_cap) {
                        j = a->head;
                        if (!j->parent) {
                            j->is_sink = 0;
                            j->parent = a->sister;
                            j->TS = i->TS;
                            j->DIST = i->DIST + 1;
                            set_active(j);
                        } else if (j->is_sink) break;
                        else if (j->TS <= i->TS &&
                                 j->DIST > i->DIST) {
                            /* heuristic - trying to make the distance from j to the source shorter */
                            j->parent = a->sister;
                            j->TS = i->TS;
                            j->DIST = i->DIST + 1;
                        }
                    }
            } else a = NULL;

            TIME++;

            if (a) {
                i->next = i; /* set active flag */
                current_node = i;

                /* augmentation */
                augment(a);
                /* augmentation end */

                /* adoption */
                while (np = orphan_first) {
                    np_next = np->next;
                    np->next = NULL;

                    while (np = orphan_first) {
                        orphan_first = np->next;
                        i = np->ptr;
                        nodeptr_block->Delete(np);
                        if (!orphan_first) orphan_last = NULL;
                        if (i->is_sink) process_sink_orphan(i);
                        else process_source_orphan(i);
                    }

                    orphan_first = np_next;
                }
                /* adoption end */
            } else current_node = NULL;
        }

        delete nodeptr_block;

        return flow;
    }

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

/* internal variables and functions */

    struct arc_st;

/* node structure */
    typedef struct node_st {
        arc_st *first;    /* first outcoming arc */

        arc_st *parent;  /* node's parent */
        node_st *next;    /* pointer to the next active node
									   (or to itself if it is the last node in the list) */
        int TS;      /* timestamp showing when DIST was computed */
        int DIST;    /* distance to the terminal */
        short is_sink;  /* flag showing whether the node is in the source or in the sink tree */

        captype tr_cap;    /* if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
									   otherwise         -tr_cap is residual capacity of the arc node->SINK */
#ifndef _METRIC_DISTANCE_
        short conflict_time;
#endif
    } node;

/* arc structure */
    typedef struct arc_st           // arc pq
    {
        node_st *head;    /* node q, i.e. node the arc points to */
        arc_st *next;    /* next arc with the same originating node */
        arc_st *sister;  /* arc qp, i.e. reverse arc */

        captype r_cap;    /* residual capacity */

        Real cap;        // cap_{pq}
    } arc;

/* 'pointer to node' structure */
    typedef struct nodeptr_st {
        node_st *ptr;
        nodeptr_st *next;
    } nodeptr;

    Block<node> *node_block;
    Block<arc> *arc_block;
    DBlock<nodeptr> *nodeptr_block;

    void (*error_function)(char *);  /* this function is called if a error occurs,
										   with a corresponding error message
										   (or exit(1) is called if it's NULL) */

    flowtype flow;        /* total flow */

/* Constructor. Optional argument is the pointer to the
   function which will be called if an error occurs;
   an error message is passed to this function. If this
   argument is omitted, exit(1) will be called. */
    Graph(node *nodes, arc *arcs, int num_nodes, void (*err_function)(char *) = NULL) {
        error_function = err_function;
        _nodes = nodes;
        _arcs = arcs;
        _num_nodes = num_nodes;
        flow = 0;
    }

    void reset_flow(void) {
        flow = 0;
    }

    arc *_arcs;
    node *_nodes;
    int _num_nodes;

//private:

/***********************************************************************/

    node *queue_first[2], *queue_last[2];  /* list of active nodes */
    nodeptr *orphan_first, *orphan_last;    /* list of pointers to orphans */
    int TIME;                /* monotonically increasing global counter */

/***********************************************************************/

/*
Functions for processing active list.
i->next points to the next node in the list
(or to i, if i is the last node in the list).
If i->next is NULL iff i is not in the list.

There are two queues. Active nodes are added
to the end of the second queue and read from
the front of the first queue. If the first queue
is empty, it is replaced by the second queue
(and the second queue becomes empty).
*/
    inline void set_active(node *i) {
        if (!i->next) {
            /* it's not in the list yet */
            if (queue_last[1]) queue_last[1]->next = i;
            else queue_first[1] = i;
            queue_last[1] = i;
            i->next = i;
        }
    }

    node *next_active() {
        node *i;

        while (1) {
            if (!(i = queue_first[0])) {
                queue_first[0] = i = queue_first[1];
                queue_last[0] = queue_last[1];
                queue_first[1] = NULL;
                queue_last[1] = NULL;
                if (!i) return NULL;
            }

            /* remove it from the active list */
            if (i->next == i) queue_first[0] = queue_last[0] = NULL;
            else queue_first[0] = i->next;
            i->next = NULL;

            /* a node in the list is active iff it has a parent */
            if (i->parent) return i;
        }
    }

    void maxflow_init() {
        node *i;

        queue_first[0] = queue_last[0] = NULL;
        queue_first[1] = queue_last[1] = NULL;
        orphan_first = NULL;

        int k;
        for (k = 0, i = _nodes; k < _num_nodes; k++, i++) {
            i->next = NULL;
            i->TS = 0;
            if (i->tr_cap > 0) {
                /* i is connected to the source */
                i->is_sink = 0;
                i->parent = TERMINAL;
                set_active(i);
                i->TS = 0;
                i->DIST = 1;
            } else if (i->tr_cap < 0) {
                /* i is connected to the sink */
                i->is_sink = 1;
                i->parent = TERMINAL;
                i->TS = 0;
                i->DIST = 1;
            } else {
                i->parent = NULL;
            }
        }
        TIME = 0;
    }

    void augment(arc *middle_arc) {
        node *i;
        arc *a;
        captype bottleneck;
        nodeptr *np;


/* 1. Finding bottleneck capacity */
/* 1a - the source tree */
        bottleneck = middle_arc->r_cap;
        for (i = middle_arc->sister->head;; i = a->head) {
            a = i->parent;
            if (a == TERMINAL) break;
            if (bottleneck > a->sister->r_cap) bottleneck = a->sister->r_cap;
        }
        if (bottleneck > i->tr_cap) bottleneck = i->tr_cap;
/* 1b - the sink tree */
        for (i = middle_arc->head;; i = a->head) {
            a = i->parent;
            if (a == TERMINAL) break;
            if (bottleneck > a->r_cap) bottleneck = a->r_cap;
        }
        if (bottleneck > -i->tr_cap) bottleneck = -i->tr_cap;


/* 2. Augmenting */
/* 2a - the source tree */
        middle_arc->sister->r_cap += bottleneck;
        middle_arc->r_cap -= bottleneck;
        for (i = middle_arc->sister->head;; i = a->head) {
            a = i->parent;
            if (a == TERMINAL) break;
            a->r_cap += bottleneck;
            a->sister->r_cap -= bottleneck;
            if (!a->sister->r_cap) {
                /* add i to the adoption list */
                i->parent = ORPHAN;
                np = nodeptr_block->New();
                np->ptr = i;
                np->next = orphan_first;
                orphan_first = np;
            }
        }
        i->tr_cap -= bottleneck;
        if (!i->tr_cap) {
            /* add i to the adoption list */
            i->parent = ORPHAN;
            np = nodeptr_block->New();
            np->ptr = i;
            np->next = orphan_first;
            orphan_first = np;
        }
/* 2b - the sink tree */
        for (i = middle_arc->head;; i = a->head) {
            a = i->parent;
            if (a == TERMINAL) break;
            a->sister->r_cap += bottleneck;
            a->r_cap -= bottleneck;
            if (!a->r_cap) {
                /* add i to the adoption list */
                i->parent = ORPHAN;
                np = nodeptr_block->New();
                np->ptr = i;
                np->next = orphan_first;
                orphan_first = np;
            }
        }
        i->tr_cap += bottleneck;
        if (!i->tr_cap) {
            i->parent = NULL;
        }

        flow += bottleneck;
    }

    void process_source_orphan(node *i) {
        node *j;
        arc *a0, *a0_min = NULL, *a;
        nodeptr *np;
        int d, d_min = INFINITE_D;

/* trying to find a new parent */
        for (a0 = i->first; a0; a0 = a0->next)
            if (a0->sister->r_cap) {
                j = a0->head;
                if (!j->is_sink && (a = j->parent)) {
                    /* checking the origin of j */
                    d = 0;
                    while (1) {
                        if (j->TS == TIME) {
                            d += j->DIST;
                            break;
                        }
                        a = j->parent;
                        d++;
                        if (a == TERMINAL) {
                            j->TS = TIME;
                            j->DIST = 1;
                            break;
                        }
                        if (a == ORPHAN) {
                            d = INFINITE_D;
                            break;
                        }
                        j = a->head;
                    }
                    if (d < INFINITE_D) /* j originates from the source - done */
                    {
                        if (d < d_min) {
                            a0_min = a0;
                            d_min = d;
                        }
                        /* set marks along the path */
                        for (j = a0->head; j->TS != TIME; j = j->parent->head) {
                            j->TS = TIME;
                            j->DIST = d--;
                        }
                    }
                }
            }

        if (i->parent = a0_min) {
            i->TS = TIME;
            i->DIST = d_min + 1;
        } else {
            /* no parent is found */
            i->TS = 0;

            /* process neighbors */
            for (a0 = i->first; a0; a0 = a0->next) {
                j = a0->head;
                if (!j->is_sink && (a = j->parent)) {
                    if (a0->sister->r_cap) set_active(j);
                    if (a != TERMINAL && a != ORPHAN && a->head == i) {
                        /* add j to the adoption list */
                        j->parent = ORPHAN;
                        np = nodeptr_block->New();
                        np->ptr = j;
                        if (orphan_last) orphan_last->next = np;
                        else orphan_first = np;
                        orphan_last = np;
                        np->next = NULL;
                    }
                }
            }
        }
    }

    void process_sink_orphan(node *i) {
        node *j;
        arc *a0, *a0_min = NULL, *a;
        nodeptr *np;
        int d, d_min = INFINITE_D;

/* trying to find a new parent */
        for (a0 = i->first; a0; a0 = a0->next)
            if (a0->r_cap) {
                j = a0->head;
                if (j->is_sink && (a = j->parent)) {
                    /* checking the origin of j */
                    d = 0;
                    while (1) {
                        if (j->TS == TIME) {
                            d += j->DIST;
                            break;
                        }
                        a = j->parent;
                        d++;
                        if (a == TERMINAL) {
                            j->TS = TIME;
                            j->DIST = 1;
                            break;
                        }
                        if (a == ORPHAN) {
                            d = INFINITE_D;
                            break;
                        }
                        j = a->head;
                    }
                    if (d < INFINITE_D) /* j originates from the sink - done */
                    {
                        if (d < d_min) {
                            a0_min = a0;
                            d_min = d;
                        }
                        /* set marks along the path */
                        for (j = a0->head; j->TS != TIME; j = j->parent->head) {
                            j->TS = TIME;
                            j->DIST = d--;
                        }
                    }
                }
            }

        if (i->parent = a0_min) {
            i->TS = TIME;
            i->DIST = d_min + 1;
        } else {
            /* no parent is found */
            i->TS = 0;

            /* process neighbors */
            for (a0 = i->first; a0; a0 = a0->next) {
                j = a0->head;
                if (j->is_sink && (a = j->parent)) {
                    if (a0->r_cap) set_active(j);
                    if (a != TERMINAL && a != ORPHAN && a->head == i) {
                        /* add j to the adoption list */
                        j->parent = ORPHAN;
                        np = nodeptr_block->New();
                        np->ptr = j;
                        if (orphan_last) orphan_last->next = np;
                        else orphan_first = np;
                        orphan_last = np;
                        np->next = NULL;
                    }
                }
            }
        }
    }
};

#endif
