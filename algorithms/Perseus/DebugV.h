/*
 * DebugV.h
 * Collection of boolean flags for debugging purposes only
 */

#ifndef DEBUGV_H_
#define DEBUGV_H_

# define SHUTUP false // if this is true then there is no debugging ouput

# define TIMEINFO (false && !SHUTUP) // which function is taking time?

# define MOVTALK (false && !SHUTUP)

# define FFTALK (false && !SHUTUP) // free face collapse talk
# define MAKEBPS (false && !SHUTUP) // adds cin.get()s strategically to interrupt flow
# define FLOWTALK (true && !SHUTUP) // program flow: which function am i in?

# define DCTTALK (false && !SHUTUP) // dense cubical toplex construction talk
# define DCFACETALK (false && !SHUTUP) // face construction for toplex
# define CDTALK (false && !SHUTUP) // dimension capping

# define BARYTALK (false && !SHUTUP) // barycentric subdivision talk!

# define UCTALK (false && !SHUTUP) // marking uncritical talk

# define DISTMATTALK (false && !SHUTUP) // talk about building rips complex from distance info
# define RIPSTALK (false && !SHUTUP) // talk about constructing rips complex
# define RIPSNBRTALK (false && !SHUTUP) // talk about constructing neighborhood matrix for rips complex
# define RIPSINTTALK (false && !SHUTUP)
# define RIPSWITTALK (false && !SHUTUP) // talk about witnessing rips complex
# define BIGRIPSTALK (false && !SHUTUP) // only top level talk about rips complex
# define INTRIPSTALK (false && !SHUTUP)
# define DIMSIMPSTALK (false && !SHUTUP)

# define BIGMORSTALK (false && !SHUTUP) // top level morse reduction
# define CHANDTALK (false && !SHUTUP)

# define CTALK (false && !SHUTUP) // talk about complexes
# define EGTALK (false && !SHUTUP) // talk about making example cell complexes
# define CMFTALK (false && !SHUTUP) // talk about iterating over coreduce

# define CORETALK (false && !SHUTUP) // coreduction algorithm debug madness
# define MTALKBIG (false && !SHUTUP) // important morse theory talk

# define KGTALK (false && !SHUTUP) // talk about storing king chains as generator completions of aces
# define GPTALK (false && !SHUTUP) // talk about gradient path computations
# define AKQTALK (false && !SHUTUP) // other coreduction talk
# define GPMEMTALK (false && !SHUTUP) // memory management for gpath cache
# define DCTSPEEDTALK (false && !SHUTUP) // speed of insertion from dense cubical toplex to cell complex

# define TOPTALK (false && !SHUTUP) // toplex talk
# define SIMTOPTALK (false && !SHUTUP) // to build simpl compex from top cells
# define SIMFACETALK (false && !SHUTUP) // to build faces of simplicial toplex
# define CUBTOPTALK (false && !SHUTUP) // to build cubical complex out of dim + anchors
# define CUBFACETALK (false && !SHUTUP) // faces of cubes!

# define PERSTALK (false && !SHUTUP) // persistence talk
# define RPRTALK (false && !SHUTUP) // persistence talk

#endif /* DEBUGV_H_ */
