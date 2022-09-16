# expand.ellipsis returns the longest possible match to the sequence.
# In the case of ambiguity, the functions stop with an error and prints at least to indistinguishable sequences that matches to the sequence
samples=list(
sam1 = 'x1,...,x5',
sam2 = 'x11,...,x51',
sam3 = 'y1a,y2b,...,y5e',
sam4 = 'x1,x3,...,x9',
sam5 = '(x11,x12,x13),(x12,x13,x14),...,(x81,x82,x83)',  #in defining vector-valued distributions in the family: block
sam6 = '(x11,x12,x13) for y1,...,(x81,x82,x83) for y8',  #in defining latent variables within a family:block
sam7 = 'x1 for (z1,y1),...,x8 for (z8,y8)',  #as in zero-inflated distributions
#By naming latent variables in a family: block, you can use both latent variable and the original variables in the model
#You can use any keyword (=>|->|as|=~|=|~|as) instead of 'for' to define latent variables
#You can use different templates for naming latent variables in a family: block
# Different templates are defined for your convenience. You can pick any template you can better to remember
sam8 = 'x1 for y1,x2 for y2,...,x8 for y8', 
sam9 = '(x11,x12 for y21,y22),(x21,x22 for y11,y12),...,(x81,x82 for y81,y82)', 
sam10 = '(x11,x12) for (y21,y22),(x21,x22) for (y11,y12),...,(x81,x82) for (y81,y82)', 
#Ellipsis can be used in regression blocks and measrement blocks as well. 
#You can premultiplier mechanism to assign labels or transform variables
sam11 = 'beta1*scale(x1)+...+beta9*scale(x9)'
)

for(i in seq_along(samples)){

  expa=expand.ellipsis(samples[[i]],details = TRUE)
  cat("\n\nSequence ",i,":\n",samples[[i]],"\nexpansion: (Length=",attr(expa,"len"),")\n")
  cat(expa)
}



