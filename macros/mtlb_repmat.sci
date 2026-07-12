function B = mtlb_repmat(A, varargin)
// Runtime emulation of Matlab's repmat(), forwarded to Scilab's native repmat().
//
// arfit was ported through Scilab's old m2sci Matlab->Scilab converter, which emits
// calls to this exact name whenever it cannot statically resolve an argument's type
// at conversion time (see modules/m2sci/macros/sci_files/sci_repmat.sci in Scilab
// core, whose own comment reads "Emulation function: mtlb_repmat()"). That file is
// only the CONVERTER rule, though -- Scilab core never shipped the matching runtime
// macro, unlike its mtlb_sum/mtlb_triu/mtlb_mean siblings in
// modules/m2sci/macros/compat_functions/. arsim.sci and arres.sci both call
// mtlb_repmat(...) on this missing shim, which raised an uncaught "Undefined
// variable: mtlb_repmat" -- and under a non-interactive scilab-adv-cli -f run (as
// used by the toolbox verification harness), an uncaught top-level error like that
// leaves the process blocked reading from a non-TTY stdin instead of exiting,
// observed externally as a hang.
//
// Scilab's native repmat() already accepts the exact calling conventions used at
// every call site here (repmat(A,[d1,d2,...]) and repmat(A,m,n)), so this is a
// straight pass-through.
    B = repmat(A, varargin(:));
endfunction
