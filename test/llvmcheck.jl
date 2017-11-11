# Functions that perform checks on the llvm rendering of methods.

"""
    creates_GC_root(f, types)

Check if a mehtod creates any garbage collection (GC) roots.
This can be used to verify functions that are expected to be called in tight
loops and should not create garbage collected objects.
"""
# Honestly, I have no idea how to write this, but "jl_get_ptls_states"
# seems to be a good shibboleth.
creates_GC_root(f, types) = contains(sprint(code_llvm, f, types), "jl_get_ptls_states")
