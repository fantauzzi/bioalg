# Code adapted from Wikipedia

####################################################
############  START MODIFICATIONS HERE  ############
####################################################

states = ('Fair', 'Biased')
 
observations = ('tails', 'heads', 'heads', 'heads', 'heads')

start_probability = {'Fair': 0.50, 'Biased': 0.50}
 
transition_probability = {
   'Fair' : {'Fair': 0.90, 'Biased': 0.10},
   'Biased' : {'Fair': 0.10, 'Biased': 0.90}
   }
 
emission_probability = {
   'Fair' : {'heads': 0.50, 'tails': 0.50},
   'Biased' : {'heads': 0.75, 'tails': 0.25}
   }
   
###################################################
############  END MODIFICATIONS HERE  #############
###################################################  

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]
 
    # Run Viterbi for t > 0
    for t in range(1, len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = max((V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states)
            V[t][y] = prob
            newpath[y] = path[state] + [y]
 
        # Don't need to remember the old paths
        path = newpath
    n = 0           # if only one element is observed max is sought in the initialization values
    if len(obs) != 1:
        n = t
    print_dptable(V)
    (prob, state) = max((V[n][y], y) for y in states)
    return (prob, path[state])
 
# Don't study this, it just prints a table of the steps.
def print_dptable(V):
    s = "    " + " ".join(("%7d" % i) for i in range(len(V))) + "\n"
    for y in V[0]:
        s += "%.5s: " % y
        s += " ".join("%.7s" % ("%f" % v[y]) for v in V)
        s += "\n"
    print(s)

def example():
    return viterbi(observations,
                   states,
                   start_probability,
                   transition_probability,
                   emission_probability)
res = example()
print(res)

for item in res[1]:
    s='E' if item =='exon' else 'I' if item=='intron' else 'A' if item in ('acceptor1', 'acceptor2') else 'D' if item in ('donor1', 'donor2') else '?'
    print(s, end='')
print()
