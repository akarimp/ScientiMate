function [d] = similaritymeasure(x, y, CalcMethod, dispout)
%{
.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2020-01-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

similaritymeasure
=================

.. code:: MATLAB

    [d] = similaritymeasure(x, y, CalcMethod, dispout)

Description
-----------

Measure similarity between two arrays

Inputs
------

x
    First array, its similarity is measured against y array
y
    Second array, its similarity is measured against x array
CalcMethod='euclidean'
    | Similarity  calculation method
    | 'euclidean': Euclidean distance
    | 'manhattan': Manhattan distance
    | 'minkowski': Minkowski distance (power=3)
    | 'cosine': Cosine distance
    | 'pearson': Pearson's correlation coefficient
    | 'spearman': spearman's correlation coefficient
    | 'norm': Absolute difference of norm
    | 'covariance': Covariance
    | 'inv_covariance': Euclidean distance of inverse covariance
    | 'histogram': Mean of absolute difference of histogram
    | 't-test': Two-sample t-test statistic
dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

d
    Arrays similarity measure

Examples
--------

.. code:: MATLAB

    x=[0,2,4,6];
    y=[2,3,5,7];
    [d]=similaritymeasure(x,y,'euclidean','yes');

    x=[1,2,3,4,5,6,7,8,9,10];
    y=[1.1,1.98,3.3,4.2,4.8,5.95,7.5,7.7,8.99,10.5];
    [d]=similaritymeasure(x,y,'pearson','yes');

References
----------
Kianimajd, A., Ruano, M. G., Carvalho, P., Henriques, J., Rocha, T., Paredes, S., & Ruano, A. E. (2017).
Comparison of different methods of measuring similarity in physiologic time series.
IFAC-PapersOnLine, 50(1), 11005-11010.

* https://en.wikipedia.org/wiki/Similarity_measure
* https://en.wikipedia.org/wiki/Goodness_of_fit
* https://dataaspirant.com/2015/04/11/five-most-popular-similarity-measures-implementation-in-python/
* https://towardsdatascience.com/similarity-measures-e3dbd4e58660
* https://www.mathworks.com/matlabcentral/answers/377944-how-to-calculate-a-percentage-of-similarity-between-two-arrays
* https://en.wikipedia.org/wiki/Template_matching
* https://www.mathworks.com/help/images/ref/normxcorr2.html
* https://www.mathworks.com/help/stats/hypothesis-tests-1.html

.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2020 Arash Karimpour
..
.. http://www.arashkarimpour.com
..
.. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
.. IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
.. FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
.. AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
.. LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
.. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
.. SOFTWARE.
%}

%--------------------------------------------------------------------------
%CODE
%--------------------------------------------------------------------------
%Assign default values

switch nargin
case 2
    CalcMethod='euclidean'; dispout='no';
case 3
    dispout='no';
end

%--------------------------------------------------------------------------
%Checking if inputs are column vectors

%if isrow(x)==1
%    x=x';
%end

%if isrow(y)==1
%    y=y';
%end

%--------------------------------------------------------------------------
%Removing NaN and Inf values

%Finding NaN and Inf in x dataset
x(((isnan(x)==1) | (isinf(x)==1)))=0;

%Finding NaN and Inf in y dataset
y(((isnan(x)==1) | (isinf(x)==1)))=0;

%--------------------------------------------------------------------------
%Calculating statistical properties for x and y data

xmean=mean(x(:)); %Mean value of x data
xstd=std(x(:)); %Standard deviation of x data

ymean=mean(y(:)); %Mean value of y data
ystd=std(y(:)); %Standard deviation of y data

%--------------------------------------------------------------------------
%Calculate similarity

%Calculate similarity using Euclidean distance
if strcmp(CalcMethod,'euclidean')==1
    %Or d=sqrt(sum((x(:)-y(:)).^2));
    p=2;
    d=(sum((abs(x(:)-y(:))).^p))^(1/p);

%Calculate similarity using Manhattan distance
elseif strcmp(CalcMethod,'manhattan')==1
    %Or d=sum(abs(x(:)-y(:)));
    p=1;
    d=(sum((abs(x(:)-y(:))).^p))^(1/p);

%Calculate similarity using Minkowski distance
elseif strcmp(CalcMethod,'minkowski')==1
    p=3;
    d=(sum((abs(x(:)-y(:))).^p))^(1/p);

%Calculate similarity using Cosine distance
elseif strcmp(CalcMethod,'cosine')==1
    d=(sum(x(:).*y(:)))/(sqrt(sum(x(:).^2))*sqrt(sum(y(:).^2)));

%Pearson's correlation coefficient
elseif strcmp(CalcMethod,'pearson')==1
    %Or d=(sum((x(:)-xmean).*(y(:)-ymean)))/(sqrt((sum((x(:)-xmean).^2))*(sum((y(:)-ymean).^2)))); %Pearson's correlation coefficient
    %Or [d, p_val] = corr(x(:),y(:),'Type','Pearson'); %Pearson's correlation coefficient
    %Or d = corr(x(:),y(:)); %Pearson's correlation coefficient
    corr_matrix = corrcoef(x(:),y(:)); %Pearson's correlation coefficient
    d=corr_matrix(1,2);

%Spearman's Spearman coefficient
elseif strcmp(CalcMethod,'spearman')==1
    %Or d=spearman(x(:),y(:));
    [d, p_val] = corr(x(:),y(:),'Type','Spearman'); %Spearman's correlation coefficient

%Absolute difference of norm
%https://www.mathworks.com/matlabcentral/answers/19752-metrics-for-matrices-similarity
elseif strcmp(CalcMethod,'norm')==1
    d=abs(norm(x(:))-norm(y(:)));

%Covariance
elseif strcmp(CalcMethod,'covariance')==1
    cov_matrix=cov(x(:),y(:));
    d=cov_matrix(1,2);

%Euclidean distance of inverse covariance
elseif strcmp(CalcMethod,'inv_covariance')==1
    cov_x=cov(x);
    cov_y=cov(y);
    if (ndims(x)>1) & (ndims(x)>1)
        cov_x_inv=inv(cov_x);
        cov_y_inv=inv(cov_y);
        d=sqrt(sum((cov_x_inv(:)-cov_y_inv(:)).^2)); %Euclidean distance
    elseif (numel(cov_x)==1) & (numel(cov_y)==1)
        d=sqrt(sum((1/cov_x-1/cov_y)**2)); %Euclidean distance
    else
        disp('Cannot calculate inverse covariance, inputs should have the same dimensions')
    end

%Mean of absolute difference of histogram
elseif strcmp(CalcMethod,'histogram')==1
    d=mean(abs(hist(x(:))-hist(y(:))));

%Two-sample t-test statistic
elseif strcmp(CalcMethod,'t-test')==1
    [h,p_val,ci,stats]=ttest2(x(:),y(:),'Vartype','unequal');
    d=stats.tstat;

end

%Normalize similarity value
%d=d/max(d(:));  %Normalize to [0, 1]

%--------------------------------------------------------------------------
%Displaying results

if strcmp(dispout,'yes')==1

    val=[d];
    name={'Similarity measure'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end

    subplot(1,2,1)
    x_y_data=[x(:),y(:)];
    hist(x_y_data)
    xlabel('Bins')
    ylabel('Frequency')
    legend('x','y')


    subplot(1,2,2)
    scatter(x,y)
    hold on
    plot([min(x(:,1));max(x(:,1))],[min(x(:,1));max(x(:,1))])
    xlabel('x')
    ylabel('y')
    legend('Data','1:1 Line')

end


%--------------------------------------------------------------------------
