cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
coder.varsize('A');
coder.varsize('phi');
codegen -config cfg eigenValueSphere3 -args {A,phi,1,1}