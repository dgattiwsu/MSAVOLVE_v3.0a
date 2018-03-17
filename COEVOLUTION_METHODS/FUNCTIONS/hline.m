% Copyright (c) 2001, Brandon Kuczenski
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

function lineHandles = hline(y, linePropertiesList, label, textPosition, textPropertiesList, axesHandle, isLineColorSynchronizedWithText)

% Draw horizontal lines (specified by a vector of y)
%
% Except the first argument 'y', everything else are optional. 
% To preserve input argument order, enter [] for unused input arguments
%
% label:        simple strings applies to all lines
%               cell strings applies to each line
% textPosition: [x, y] bewteen 0 and 1 relative to the 'lower' end
%
% Properties lists are unpacked as tail input arguments to text/plot():
% linePropertiesList: a cell that specify properties in plot()
%                     it can also be used to specify line types in plot()
%                     [simple string input will be converted to cellstr]
% textPropertiesList: a cell that specify properties in text()
%
% isLineColorSynchronizedWithText: by default it's true. 
%                                  Effects might override textPropertiesList
%
% Example:
%
%   figure, t = 0:0.01:2*pi; X=sin(t); Y=cos(t); plot(t', [X', Y'])
%
%   Simple usage: default plot() settings, no label, current axes
%       hline([0.5 -0.5]);
%
%   Simple line control: 
%       hline([0.5 -0.5], 'r:');
%
%   Advanced line control: 
%       hline([0.5 -0.5], {'r:', 'LineWidth', 5});
%
%   Display simple labels: default text() settings, sync colors
%       hline([0.5 -0.5], ':', {'pi/2', 'pi'});
%
%   Control text position (move label to right) and TeX display
%       hline([0.5 -0.5], ':', {'\pi/2', '\pi'}, [1 0], {'Interpreter', 'tex'});

% Author:      Hoi Wong (wonghoi.ee@gmail.com)
% Date:        02/14/2008
%
% Updated:     03/11/2011 (Synchronized new features with vline())
%
% Acknowledgement:  It was based on hline() written by Brandon Kuczenski.

    if( ~isvector(y) )
        error('y must be a vector');
    else
        y=y(:);     % Make sure it's column vector
    end

    if( ~exist('label', 'var') )                label = [];                 end
    if( ~exist('textPosition', 'var') )         textPosition = [];          end
    if( ~exist('linePropertiesList', 'var') )   linePropertiesList = {};    end
    if( ~exist('textPropertiesList', 'var') )   textPropertiesList = {};    end
    if( ~exist('axesHandle', 'var') )           axesHandle = gca;           end
    if( ~exist('isLineColorSynchronizedWithText', 'var') )    
        isLineColorSynchronizedWithText = true;    
    end
    
    if( isempty(axesHandle) )           axesHandle = gca;           end
    if( isempty(textPropertiesList) )   textPropertiesList = {};    end
    if( isempty(linePropertiesList) )   linePropertiesList = {};    end  
    
    if( ~ishandle(axesHandle) )
        error('axesHandle specified is invalid (either empty or real graphics handle please)');
    end
    if( ischar(linePropertiesList) )    linePropertiesList = {linePropertiesList};    end    
    
    if( isempty(textPosition) )
        textPositionX = 0.02;
        textPositionY = 0.02;
        % horizontal lines doesn't need text rotation like vertical lines
    elseif( isscalar(textPosition) )
        textPositionX = textPosition;
        textPositionY = 0.02; 
    elseif( length(textPosition)~=2 )
        error('Invalid textPosition');
    else
        textPositionX = textPosition(1);
        textPositionY = textPosition(2);
    end
    clear textPosition;     % Avoid typos for similarly named variables
               
    holdState = ishold(axesHandle);
    hold(axesHandle, 'on');
    
    xLimits=get(axesHandle,'xlim');             % Row vector    
    Xlimits=repmat(xLimits', 1, length(y));
    
    % Example: for horizontal lines
    % X = [2 2        Y = [3 4
    %      5 5];           3 4];
            
    lineHandles = plot(axesHandle, Xlimits, [y';y'], linePropertiesList{:});
    
    if( ~isempty(label) )
        yLimits = get(axesHandle,'ylim');
        yLowerLimit = yLimits(2);
        yUpperLimit = yLimits(1);        
        yRange      = yUpperLimit - yLowerLimit;
        yPosition   = y-textPositionY*yRange;
        
        xUpperLimit = xLimits(2);
        xLowerLimit = xLimits(1);
        xRange      = xUpperLimit - xLowerLimit;
        xPosition   = xLowerLimit+textPositionX*xRange;
        Xposition   = repmat(xPosition, length(y), 1);
        
        % textHandle is a vector correspond to the line
        textHandles = text(Xposition, yPosition, label, 'Parent', axesHandle, textPropertiesList{:});
        
        % Set the text colors to be identical to line colors
        if( isLineColorSynchronizedWithText )
            arrayfun(@(k) set( textHandles(k), 'color', get(lineHandles(k), 'color') ), 1:length(y));
        end
                 
    end
    
    if( holdState==false )
        hold(axesHandle, 'off');
    end
    % this last part is so that it doesn't show up on legends
    set(lineHandles,'tag','hline','handlevisibility','off') 

