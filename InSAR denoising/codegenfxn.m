cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
codegen -config cfg InSAR_DG -args {g,p,q,alpha,tol,xtol,dt}