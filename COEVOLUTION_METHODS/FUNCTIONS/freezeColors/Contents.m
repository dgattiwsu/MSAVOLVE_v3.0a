% Copyright (c) 2017, John Iversen
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Functions enabling use of multiple colormaps per figure.
%
% version 2.3, 3/2007
%
%   freezeColors    - Lock colors of plot, enabling multiple colormaps per figure.
%   unfreezeColors  - Restore colors of a plot to original indexed color.
%   
% demo/html
%   freezeColors_pub.html   - Published demonstration.
%
% test
%   test_main       - Test function.

% AUTHOR
% John Iversen, 2005-10
% john_iversen@post.harvard.edu
% 
% Free for any use, so long as AUTHOR information remains with code.
%
% HISTORY
% 
% JRI 6/2005  (v1)
% JRI 9/2006  (v2.1) now operates on all objects with CData (not just images as before)
% JRI 11/2006 handle NaNs in images/surfaces (was not posted on file exchange, included in v2.3)
% JRI 3/2007  (v2.3) Fix bug in handling of caxis that made colorbar scales on frozen plots incorrect. 
% JRI 4/2007  Add 'published' html documentation.
% JRI 9/2010  Changed documentation for colorbars
