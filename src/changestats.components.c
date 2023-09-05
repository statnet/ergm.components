/*
######################################################################
#
# changestats.components.c
#
# copyright (c) 2022, Carter T. Butts <buttsc@uci.edu>
# Last Modified 08/02/23
# Licensed under the GNU General Public License v3 or later
#
# Portions based on the statnet ergm.userterms package (Hunter et al.,
# 2013, JSS).
#
# Part of the R/ergm.components package
#
# This file contains routines related to the changescores for
# component and connectivity-related ergm terms.
#
######################################################################
*/

/*HEADERS--------------------------------------------------------------------*/

#include "changestats.components.h"


/*UTILITY FUNCTIONS----------------------------------------------------------*/

/*
Local (recursive) depth-first search function used by localBridgeCnt, to find
bridging edges involving the component(s) containing head and tail.  The search
is started at v; root is the current search root, thpresent indicates whether
the {tail,head} edge should be treated as present vs. absent, and iter,
visited, miniter, and visiter are pointers to memory structures containing the
DFS state.  bc is a pointer to the bridge count.  *bc is updated as bridges
are found.

This function will work with undirected networks, or with directed networks
when semi-bridging dyads are to be counted.  In the latter case, we actually
count semi-bridging dyads when dyadic==TRUE, and otherwise we count weak
cutedges; the difference is that a weak cutedge must both be part of a semi-
bridge and also have no reciprocating edge.

Note that, since the changescore is zero in the directed, dyadic case where
(head,tail) is in the graph, upstream code will not/should not call this routine
in that case.  Thus, we are free to assume here that we are only being called
in the directed/dyadic case where there is no (head,tail) edge.  If using this
code for something else, be aware that it may employ that assumption.
*/
void lBCDFS(Vertex v, Vertex root, Network *nwp, Vertex tail, Vertex head, int thpresent, int dyadic, Vertex *iter, char *visited, Vertex *miniter, Vertex *viter, double *bc){
  Edge e;
  Vertex alter;
  int thcheck=0,recip=0,thumatch=0,thdmatch=0,htdmatch=0;
  
  //Rprintf("Pushing %ld (thpres=%d, ht={%ld,%ld}, parent=%ld)\n",v,thpresent,tail,head,root);
  visited[v-1]++;
  miniter[v-1]=*iter;
  viter[v-1]=(*iter)++;

  STEP_THROUGH_INEDGES(v,e,alter){
    thumatch=DYADMATCH(tail,head,alter,v);        //{tail,head} == {alter,v}
    thdmatch=DDYADMATCH(tail,head,alter,v);       //(tail,head) == (alter,v)
    //Rprintf("  Cand edge %ld=>%ld\n",v,alter);
    //Process if alter is not the node that nominated v, and:
    //  - Undirected case: ({tail,head} != {v,alter}) or thpresent 
    //  - Directed case: ((v,alter) != (head,tail)) or thpresent
    if((thpresent||(DIRECTED&&(!thdmatch))||((!DIRECTED)&&(!thumatch)))&&(alter!=root)){
      if(DYADMATCH(tail,head,v,alter)) //{v,alter} == {tail,head} (so we saw {t,h})
        thcheck++;
      //Rprintf("  %ld->%ld\n",v,alter);
      if(visited[alter-1])   /*If seen before, update minvis for v*/
        miniter[v-1]=MIN(miniter[v-1],viter[alter-1]);
      else{
        //Rprintf("    %ld is new - recursing\n",alter);
        lBCDFS(alter, v, nwp, tail, head, thpresent, dyadic, iter, visited, miniter, viter, bc);
        miniter[v-1]=MIN(miniter[v-1],miniter[alter-1]);
        if(DIRECTED&&(!dyadic)){ /*Is there a reciprocating edge?  Only care for dir cutedges*/
          recip=(DDYADMATCH(v,alter,tail,head) ? thpresent : IS_OUTEDGE(v,alter));
        }
        //Increment count if miniter[alter-1]>viter[v-1] and undirected or no reciprocator
        if((miniter[alter-1]>viter[v-1])&&((!DIRECTED)||dyadic||(!recip)))
          (*bc)++;
      }
    }
  }
  STEP_THROUGH_OUTEDGES(v,e,alter){
    thumatch=DYADMATCH(tail,head,v,alter);        //{tail,head} == {v,alter}
    thdmatch=DDYADMATCH(tail,head,v,alter);       //(tail,head) == (v,alter)
    htdmatch=DDYADMATCH(tail,head,alter,v);       //(tail,head) == (alter,v)
    //We only want to traverse each dyad once.  In the undirected case, that's not
    //an issue: every pair is only encoded in one direction.  In the directed case,
    //we must be careful.  We already considered the alter->v ties, so we now want
    //to rule out cases where (1) there is an alter->v normal edge, or (2) (alter,v)
    //is actually (tail,head) and thpresent==TRUE.  In either of those cases, we
    //already went down this path....
    if((!DIRECTED)||(DIRECTED&&((htdmatch&&(!thpresent))||((!htdmatch)&&(!IS_OUTEDGE(alter,v)))))){       /*Ensure we only traverse the dyad once*/
      //Rprintf("  Cand edge %ld=>%ld\n",v,alter);
      //Process if alter is not the node that nominated v, and:
      //  - Undirected case: ({tail,head} != {v,alter}) or thpresent 
      //  - Directed case: ((v,alter) != (tail,head)) or thpresent
      if((thpresent||(DIRECTED&&(!thdmatch))||((!DIRECTED)&&(!thdmatch)))&&(alter!=root)){
        if(thdmatch)    //We saw the {tail,head} dyad
          thcheck++;
        //Rprintf("  %ld->%ld\n",v,alter);
        if(visited[alter-1])   /*If seen before, update minvis for v*/
          miniter[v-1]=MIN(miniter[v-1],viter[alter-1]);
        else{
          //Rprintf("    %ld is new - recursing\n",alter);
          lBCDFS(alter, v, nwp, tail, head, thpresent, dyadic, iter, visited, miniter, viter, bc);
          miniter[v-1]=MIN(miniter[v-1],miniter[alter-1]);
          if(DIRECTED&&(!dyadic)){ /*Is there a reciprocating edge?  Only care for dir cutedges*/
            recip=(htdmatch ? thpresent : IS_OUTEDGE(alter,v));
          }
          //Increment count if miniter[alter-1]>viter[v-1] and undirected or no reciprocator
          if((miniter[alter-1]>viter[v-1])&&((!DIRECTED)||dyadic||(!recip)))
            (*bc)++;
        }
      }
    }
  }
  if(thpresent&&(!thcheck)){  /*If {tail,head} ought to be here and we didn't see it, deal w/it*/
    if((v==tail)&&(head!=root)){
      alter=head;
      thcheck=1;
    }else if((v==head)&&(tail!=root)){
      alter=tail;
      thcheck=1;
    }else
      thcheck=0;
    if(thcheck){
      //Rprintf("  %ld->%ld\n",v,alter);
      if(visited[alter-1])   /*If seen before, update minvis for v*/
        miniter[v-1]=MIN(miniter[v-1],viter[alter-1]);
      else{
        //Rprintf("    %ld is new - recursing\n",alter);
        lBCDFS(alter, v, nwp, tail, head, thpresent, dyadic, iter, visited, miniter, viter, bc);
        miniter[v-1]=MIN(miniter[v-1],miniter[alter-1]);
        if(miniter[alter-1]>viter[v-1]) //If we got here, no reciprocating edge....
          (*bc)++;
      }
    }
  }
  //Rprintf("Leaving %ld\n",v);
}


/*
Function to return number of bridges involving the component(s) containing head
and tail.  If thpresent, we assume that {tail,head} is present, else that it is
absent.

This function will work with undirected networks, or with directed networks
when semi-bridging dyads are to be counted.  In the latter case, we actually
count semi-bridging dyads when dyadic==TRUE, and otherwise we count weak
cutedges; the difference is that a weak cutedge must both be part of a semi-
bridge and also have no reciprocating edge.

FWIW, we will never call this function in the dyadic case when (head,tail) is in
the graph, since in that case the changescore would be zero.  So don't reuse this
code in an application that doesn't depend on that assumption, since it may have
been woven in, in some way that will surprise you...
*/
double localBridgeCnt(Network *nwp, Vertex tail, Vertex head, int thpresent, int dyadic){
  double bc=0.0;
  char *visited;
  Vertex *miniter,*viter,*iter;

  /*Initialize*/
  visited=Calloc(N_NODES,char);
  miniter=Calloc(N_NODES,Vertex);
  viter=Calloc(N_NODES,Vertex);
  iter=Calloc(1,Vertex);
  *iter=1;

  /*If this is directed and dyadic, the changescore is always zero if there is an incoming
    (head,tail) edge; thus only compute if we aren't in that situation.  But honestly, we
    shouldn't be called in that situation anyway, because the upstream code is not supposed
    to do so.  This is a redundant check, but still a cheap one...*/
  if(!(DIRECTED&&dyadic&&IS_OUTEDGE(head,tail))){
    /*Perform DFS*/
    lBCDFS(tail, 0, nwp, tail, head, thpresent, dyadic, iter, visited, miniter, viter, &bc);
    if(!visited[head-1]){  /*Not in the same component*/
      //Rprintf("%ld not reached - running separately\n",head);
      lBCDFS(head, 0, nwp, tail, head, thpresent, dyadic, iter, visited, miniter, viter, &bc);
    }
  }

  /*Clean up and return*/
  Free(visited);
  Free(miniter);
  Free(viter);
  Free(iter);
  return bc;
}


/*CHANGESCORE FUNCTIONS------------------------------------------------------*/

/*
Dimer count statistic - undirected version.  Counts changes in the number of
components of order 2.
*/
CHANGESTAT_FN(d_dimers_udir) {
  Vertex tail,head,alter=0,v,hdeg,tdeg;
  Edge e;
  int i,thpres;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i);      /*Get endpoints*/
    head=HEAD(i);
    thpres=IS_UNDIRECTED_EDGE(tail,head);
    tdeg=IN_DEG[tail]+OUT_DEG[tail];  /*Get current degrees*/
    hdeg=IN_DEG[head]+OUT_DEG[head];
    if(thpres){                         /*Removing edge*/
      /*
      Edge removal cases:
        - t,h both have degree 1 - eliminates one dimer
        - t or h has degree 2, neighbor has degree 1 - create one dimer per endpoint
          making criteria
        - Nothing else matters
      */
      if((tdeg==1)&&(hdeg==1)){  /*Remove {t,h} dimer*/
        CHANGE_STAT[0]--;
      }else{                     /*Creative destruction*/
        /*Check to see if tail will now be in a dimer*/
        if(tdeg==2){
          /*Find other neighbor*/
          STEP_THROUGH_OUTEDGES(tail,e,v){  /*Select only non-head*/
            if(v!=head)
              alter=v;
          }
          STEP_THROUGH_INEDGES(tail,e,v){  /*Select only non-head*/
            if(v!=head)
              alter=v;
          }
          if(IN_DEG[alter]+OUT_DEG[alter]==1) /*If alter is pendant, dimer created*/
            CHANGE_STAT[0]++;
        }
        /*Check to see if head will now be in a dimer*/
        if(hdeg==2){
          /*Find other neighbor*/
          STEP_THROUGH_OUTEDGES(head,e,v){  /*Select only non-tail*/
            if(v!=tail)
              alter=v;
          }
          STEP_THROUGH_INEDGES(head,e,v){  /*Select only non-head*/
            if(v!=tail)
              alter=v;
          }
          if(IN_DEG[alter]+OUT_DEG[alter]==1) /*If alter is pendant, dimer created*/
            CHANGE_STAT[0]++;
        }
      }
    }else{                              /*Adding edge*/
      /*
      Edge addition cases:
        - t,h both have degree 0 - creates one dimer
        - t or h has degree 1, neighbor has degree 1 - eliminate one dimer per endpoint
          making criteria
        - Nothing else matters
      */
      if((tdeg==0)&&(hdeg==0)){  /*Add {t,h} dimer*/
        CHANGE_STAT[0]++;
      }else{                     /*Destructive creation*/
        /*Check to see if tail no longer in a dimer*/
        if(tdeg==1){
          /*Find other neighbor*/
          STEP_THROUGH_OUTEDGES(tail,e,v){  /*Select only non-head*/
            if(v!=head)
              alter=v;
          }
          STEP_THROUGH_INEDGES(tail,e,v){  /*Select only non-head*/
            if(v!=head)
              alter=v;
          }
          if(IN_DEG[alter]+OUT_DEG[alter]==1) /*If alter is pendant, dimer destroyed*/
            CHANGE_STAT[0]--;
        }
        /*Check to see if head no longer in a dimer*/
        if(hdeg==1){
          /*Find other neighbor*/
          STEP_THROUGH_OUTEDGES(head,e,v){  /*Select only non-tail*/
            if(v!=tail)
              alter=v;
          }
          STEP_THROUGH_INEDGES(head,e,v){  /*Select only non-head*/
            if(v!=tail)
              alter=v;
          }
          if(IN_DEG[alter]+OUT_DEG[alter]==1) /*If alter is pendant, dimer created*/
            CHANGE_STAT[0]--;
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Dimer count statistic - directed version.  Counts changes in the number of
semi-components of order 2.
*/
CHANGESTAT_FN(d_dimers_dir) {
  Vertex tail,head,v=-1,hid,tid,hod,tod;
  Edge e;
  int i;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i);      /*Get endpoints*/
    head=HEAD(i);
    tid=IN_DEG[tail];  /*Get current degrees*/
    hid=IN_DEG[head];
    tod=OUT_DEG[tail];
    hod=OUT_DEG[head];
    if(IS_OUTEDGE(tail,head)){  /*Removing edge*/
      if((tod==1)&&(hid==1)&&(tid+hod==0)) /*Decrement if removing a dimer*/
        CHANGE_STAT[0]--;
      else if(!IS_OUTEDGE(head,tail)){ /*Increment for each dimer created by removal*/
        /*Is tail now a pendant? (Not counting tail->head; there is no head->tail)*/
        if((tod+tid>1)&&(tod<3)&&(tid<2)){  /*Perhaps, perhaps, perhaps...*/
          if(tod>1){  /*Has a non-head out-alter, maybe also an incoming tie*/
            STEP_THROUGH_OUTEDGES(tail,e,v){
              if(v!=head){
                CHANGE_STAT[0] += (tid-IS_OUTEDGE(v,tail)==0) && (IN_DEG[v]==1) && (OUT_DEG[v]-IS_OUTEDGE(v,tail)==0);                
              }
            }
          }else{      /*No non-head out-alter, only an incoming tie (and not from head)*/
            STEP_THROUGH_INEDGES(tail,e,v){
              CHANGE_STAT[0] += (OUT_DEG[v]==1) && (IN_DEG[v]==0);
            }
          }
        }
        /*Is head now a pendant? (Not counting tail->head; there is no head->tail)*/
        if((hod+hid>1)&&(hod<2)&&(hid<3)){  /*Perhaps, perhaps, perhaps...*/
          if(hid>1){  /*Has a non-tail in-alter, maybe also an outgoing tie*/
            STEP_THROUGH_INEDGES(head,e,v){
              if(v!=tail){
                CHANGE_STAT[0] += (hod-IS_OUTEDGE(head,v)==0) && (OUT_DEG[v]==1) && (IN_DEG[v]-IS_OUTEDGE(head,v)==0);                
              }
            }
          }else{      /*No non-tail in-alter, only an outgoing tie (and not to tail)*/
            STEP_THROUGH_OUTEDGES(head,e,v){
              CHANGE_STAT[0] += (OUT_DEG[v]==0) && (IN_DEG[v]==1);
            }
          }
        }
      }
    }else{                     /*Adding edge*/
      if(tod+hid+tid+hod==0){           /*Increment if creating a dimer*/
        CHANGE_STAT[0]++;
      }else if(!IS_OUTEDGE(head,tail)){ /*Decrement for each dimer lost to aggregation*/
        /*Was tail in a dimer?*/
        if((tod+tid>0)&&(tid<2)&&(tod<2)){  /*Maybe it was...*/
          if(tod>0){  /*Has one out-edge, maybe an in-edge*/
            STEP_THROUGH_OUTEDGES(tail,e,v){
              CHANGE_STAT[0] -= ((tid==0)||IS_OUTEDGE(v,tail)) && (IN_DEG[v]==1) && (OUT_DEG[v]-IS_OUTEDGE(v,tail)==0);
            }
          }else{      /*Only has an in-edge - just verify that alter is pendant*/
            STEP_THROUGH_INEDGES(tail,e,v){
              CHANGE_STAT[0]-=(OUT_DEG[v]==1)&&(IN_DEG[v]==0);
            }
          }
        }
        /*Was head in a dimer?*/
        if((hod+hid>0)&&(hid<2)&&(hod<2)){  /*Maybe it was...*/
          if(hod>0){  /*Has one out-edge, maybe an in-edge*/
            STEP_THROUGH_OUTEDGES(head,e,v){
              CHANGE_STAT[0] -= ((hid==0)||IS_OUTEDGE(v,head)) && (IN_DEG[v]==1) && (OUT_DEG[v]-IS_OUTEDGE(v,head)==0);
            }
          }else{      /*Only has an in-edge - just verify that alter is pendant*/
            STEP_THROUGH_INEDGES(head,e,v){
              CHANGE_STAT[0]-=(OUT_DEG[v]==1)&&(IN_DEG[v]==0);
            }
          }
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Component count statistic - change in the number of components when the focal edge
is added/removed.

This function will work with undirected networks, or with directed networks
when weak components are to be counted.
*/
CHANGESTAT_FN(d_components) {
  Vertex tail,head,v,alter;
  char *visited=NULL;
  element *el,*tovisit,*nextvis;
  Edge e;
  int i,flag;

  /*
  Obviously, adding an edge can either leave the number of components unchanged, or
  decrease it by 1.  Likewise, removing an edge can either leave the number of 
  components unchanged, or increase it by 1.  The cases can be distinguished as
  follows (noting that we are considering only undirected graphs here):
    - If there is a (t,h) path not using {t,h}, then neither adding nor removing
      {t,h} has any effect.
    - If there is not a (t,h) path not using {t,h}, then the changescore is -/+ 1
      for adding/removing {t,h}.
  Thus, what we have to do to calculate the changescore is to look for that path.
  We can do this via BFS in linear time, terminating early if we get lucky. 
  */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i);      /*Get endpoints*/
    head=HEAD(i);
    //Rprintf("Working on {%ld,%ld}\n",tail,head);
    el=Calloc(1,element);
    el->v=tail;
    el->next=NULL;
    el->prev=NULL;
    tovisit=el;
    nextvis=el;
    flag=(DIRECTED)&&(IS_OUTEDGE(head,tail)); /*Redundant w/in dyad edge shortcut*/
    visited=Calloc(N_NODES,char);
    visited[tail-1]=1;
    while((!flag)&&(tovisit!=NULL)){
      /*Take the next element*/
      el=nextvis;
      nextvis=el->prev;     /*Update the bottom pointer*/
      if(el->prev==NULL)    /*If this was the top element, set top ptr to NULL*/
        tovisit=NULL;
      else
        el->prev->next=NULL;
      v=el->v;
      //Rprintf("\tDeqeued %ld\n",v);
      Free(el);             /*Recycle this element*/
      /*Process v*/
      STEP_THROUGH_OUTEDGES(v,e,alter){
        if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])&&(!flag)){
          if(alter==head)  /*If we've found our goal, set the flag*/
            flag=1;
          el=Calloc(1,element);  /*Add to the top of our visitation queue*/
          el->v=alter;
          el->next=tovisit;
          el->prev=NULL;
          if(tovisit!=NULL)
            tovisit->prev=el;
          else                    /*If queue was empty, update bottom ptr*/
            nextvis=el;
          tovisit=el;
          visited[alter-1]=1;
          //Rprintf("\t\tEnqueued %ld\n",alter);
        }
      }
      STEP_THROUGH_INEDGES(v,e,alter){
        if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])&&(!flag)){
          if(alter==head)  /*If we've found our goal, set the flag*/
            flag=1;
          el=Calloc(1,element);  /*Add to the top of our visitation queue*/
          el->v=alter;
          el->next=tovisit;
          el->prev=NULL;
          if(tovisit!=NULL)
            tovisit->prev=el;
          else                    /*If queue was empty, update bottom ptr*/
            nextvis=el;
          tovisit=el;
          visited[alter-1]=1;
          //Rprintf("\t\tEnqueued %ld\n",alter);
        }
      }
    }
    /*Perform changestat increment*/
    if(!flag){
      if(DIRECTED){
        CHANGE_STAT[0]+=(IS_OUTEDGE(tail,head)||IS_OUTEDGE(head,tail) ? 1.0 : -1.0);
      }else{
        CHANGE_STAT[0]+=(IS_UNDIRECTED_EDGE(tail,head) ? 1.0 : -1.0);
      }
    }
    /*Clean up*/
    //Rprintf("Cleaning up\n");
    Free(visited);
    el=tovisit;
    while(el!=NULL){
      //Rprintf("\tRemoving %ld\n",el->v);
      nextvis=el;
      el=el->next;
      Free(nextvis);
    }
    //Rprintf("\tDone cleaning\n");
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for sums of functions of component sizes, i.e.

  sum_c f(|V(c)|)

where the sum is over components.  At present, two classes of functions are
implemented:

  f(x) = x^k

and

  f(x) = log(x)^k.
  
INPUT_PARAM[0] turns logging on/off (1=log), and k=INPUT_PARAM[1].

This function will work with undirected networks, or with directed networks
when weak components are to be counted.
*/
CHANGESTAT_FN(d_compsizesum) {
  Vertex tail,head,v,alter;
  char *visited=NULL;
  element *el,*tovisit,*nextvis;
  Edge e;
  int i,flag;
  double csh=0.0,cst=0.0,csv=0.0;

  /*
  Obviously, adding an edge can either leave the number of components unchanged, or
  decrease it by 1.  Likewise, removing an edge can either leave the number of 
  components unchanged, or increase it by 1.  The cases can be distinguished as
  follows (noting that we are considering only undirected graphs here):
    - If there is a (t,h) path not using {t,h}, then neither adding nor removing
      {t,h} has any effect.
    - If there is not a (t,h) path not using {t,h}, then the changescore is +/-
      f(Sth) - f(St) - f(Sh) for adding/removing {t,h}.
  Thus, what we have to do to calculate the changescore is to look for that path.
  We can do this via BFS in linear time, terminating early if we get lucky.  If we
  do not, then we use the BFS to tell us how large one of the endpoints' components
  is, using a BFS on the other endpoint to get the other component size.  Obviously,
  if we merge the two (adding the edge) then we get a single component of size equal
  to the sum.
  */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i);      /*Get endpoints*/
    head=HEAD(i);
    //Rprintf("Working on {%ld,%ld}\n",tail,head);
    /*First step: perform a BFS from tail, either finding head (in which case the CS
    is 0) or else calculating the size of tail's component.*/
    el=Calloc(1,element);
    el->v=tail;
    el->next=NULL;
    el->prev=NULL;
    tovisit=el;
    nextvis=el;
    flag=(DIRECTED)&&(IS_OUTEDGE(head,tail)); /*Redundant w/in dyad edge shortcut*/
    visited=Calloc(N_NODES,char);
    visited[tail-1]=1;
    cst=1.0; /*Tail component size*/
    while((!flag)&&(tovisit!=NULL)){
      /*Take the next element*/
      el=nextvis;
      nextvis=el->prev;     /*Update the bottom pointer*/
      if(el->prev==NULL)    /*If this was the top element, set top ptr to NULL*/
        tovisit=NULL;
      else
        el->prev->next=NULL;
      v=el->v;
      //Rprintf("\tDeqeued %ld\n",v);
      Free(el);             /*Recycle this element*/
      /*Process v*/
      STEP_THROUGH_OUTEDGES(v,e,alter){
        if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])&&(!flag)){
          if(alter==head)  /*If we've found our goal, set the flag*/
            flag=1;
          el=Calloc(1,element);  /*Add to the top of our visitation queue*/
          el->v=alter;
          el->next=tovisit;
          el->prev=NULL;
          if(tovisit!=NULL)
            tovisit->prev=el;
          else                    /*If queue was empty, update bottom ptr*/
            nextvis=el;
          tovisit=el;
          visited[alter-1]=1;
          cst++;  /*Increment the tail component size*/
          //Rprintf("\t\tEnqueued %ld\n",alter);
        }
      }
      STEP_THROUGH_INEDGES(v,e,alter){
        if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])&&(!flag)){
          if(alter==head)  /*If we've found our goal, set the flag*/
            flag=1;
          el=Calloc(1,element);  /*Add to the top of our visitation queue*/
          el->v=alter;
          el->next=tovisit;
          el->prev=NULL;
          if(tovisit!=NULL)
            tovisit->prev=el;
          else                    /*If queue was empty, update bottom ptr*/
            nextvis=el;
          tovisit=el;
          visited[alter-1]=1;
          cst++;  /*Increment the tail component size*/
          //Rprintf("\t\tEnqueued %ld\n",alter);
        }
      }
    }
    Free(visited);  /*Clean up*/
    el=tovisit;
    while(el!=NULL){
      nextvis=el;
      el=el->next;
      Free(nextvis);
    }
    /*Second step: if we didn't find head, then this edge is a bridge, and we need to
    find the size of head's component (without the tail,head edge, that is).*/
    if(!flag){
      el=Calloc(1,element);
      el->v=head;
      el->next=NULL;
      el->prev=NULL;
      tovisit=el;
      nextvis=el;
      visited=Calloc(N_NODES,char);
      visited[head-1]=1;
      csh=1.0; /*Head component size*/
      while(tovisit!=NULL){
        /*Take the next element*/
        el=nextvis;
        nextvis=el->prev;     /*Update the bottom pointer*/
        if(el->prev==NULL)    /*If this was the top element, set top ptr to NULL*/
          tovisit=NULL;
        else
          el->prev->next=NULL;
        v=el->v;
        Free(el);             /*Recycle this element*/
        /*Process v*/
        STEP_THROUGH_OUTEDGES(v,e,alter){
          if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])){
            el=Calloc(1,element);  /*Add to the top of our visitation queue*/
            el->v=alter;
            el->next=tovisit;
            el->prev=NULL;
            if(tovisit!=NULL)
              tovisit->prev=el;
            else                    /*If queue was empty, update bottom ptr*/
              nextvis=el;
            tovisit=el;
            visited[alter-1]=1;
            csh++;  /*Increment the head component size*/
          }
        }
        STEP_THROUGH_INEDGES(v,e,alter){
          if((!DYADMATCH(tail,head,v,alter))&&(!visited[alter-1])){
            el=Calloc(1,element);  /*Add to the top of our visitation queue*/
            el->v=alter;
            el->next=tovisit;
            el->prev=NULL;
            if(tovisit!=NULL)
              tovisit->prev=el;
            else                    /*If queue was empty, update bottom ptr*/
              nextvis=el;
            tovisit=el;
            visited[alter-1]=1;
            csh++;  /*Increment the head component size*/
          }
        }
      }
      Free(visited);  /*Clean up*/
      el=tovisit;
      while(el!=NULL){
        nextvis=el;
        el=el->next;
        Free(nextvis);
      }
    }
    /*Perform changestat increment*/
    if(!flag){
      if(INPUT_PARAM[0]>0.0){  /*Log everything first?*/
        csv=log(csh+cst);
        csh=log(csh);
        cst=log(cst);
      }else
        csv=csh+cst;
      /*Base CS: f(csv) - f(cst) - f(csh)*/
      if(INPUT_PARAM[1]==0.0){             /*Avoid numerical oddities*/
        csv = -1.0;
      }else if(INPUT_PARAM[1]==1.0){       /*Avoid slow calls*/
        csv = csv-cst-csh;
      }else if(INPUT_PARAM[1]==2.0){
        csv = csv*csv-cst*cst-csh*csh;
      }else{
        csv = pow(csv,INPUT_PARAM[1])-pow(cst,INPUT_PARAM[1])-pow(csh,INPUT_PARAM[1]);
      }
      if(DIRECTED){
        CHANGE_STAT[0]+=(IS_OUTEDGE(tail,head)||IS_OUTEDGE(head,tail) ? -csv : csv);
      }else{
        CHANGE_STAT[0]+=(IS_UNDIRECTED_EDGE(tail,head) ? -csv : csv);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Change statistic for the number of bridging edges (i.e., cutedges).

This function will work with undirected networks, or with directed networks
when bridges are defined with respect to weak connectivity.  This can be
interpreted in two ways.  Then the dyadic option is set (INPUTS_PARAM[0]),
we count bridging dyads (i.e, dyads that are bridges in the semigraph).  
When the dyadic option is unset, we count bridging edges (i.e., cutedges
with respect to weak connectivity).  This means that, in the second case,
e.g., a T021D contains 2 bridges, but a T111D contains only 1 (since
neither edge within the mutual dyad can (weakly) disconnect it).  By
contrast, in the dyadic case, both have two bridges.  In the dyadic case,
note that the changescore for an i->j edge is always zero if there exists 
a j->i  edge, since in neither case does the removal of an edge remove a
weak bridge if it is part of a mutual dyad.
*/
CHANGESTAT_FN(d_bridges) {
  Vertex tail,head;
  double bc0,bc1;
  int i, thpresent, dyadic=INPUT_PARAM[0];

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i);      /*Get endpoints*/
    head=HEAD(i);
    thpresent = (DIRECTED ? IS_OUTEDGE(tail,head) : IS_UNDIRECTED_EDGE(tail,head));
    /*Note: if directed and h->t, then toggling t->h can't do anything in the
      dyadic case; thus, we detect this case and avoid spendng time on it....*/
    if((!DIRECTED)||(!dyadic)||(!IS_OUTEDGE(head,tail))){
      //Rprintf("Working on {%ld,%ld}\n",tail,head);
      //Rprintf("Edge absent case:\n");
      bc0=localBridgeCnt(nwp,tail,head,0,dyadic);  /*Bridges w/out edge*/
      //Rprintf("Edge present case:\n");
      bc1=localBridgeCnt(nwp,tail,head,1,dyadic);  /*Bridges with edge*/
      //Rprintf("{%ld,%ld} bridges w/out %.0f, with %.0f\n",tail,head,bc0,bc1);
      CHANGE_STAT[0]+=(thpresent ? bc0-bc1 : bc1-bc0);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}



