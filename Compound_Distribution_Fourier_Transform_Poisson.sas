 %let timenow=%sysfunc(time(), time.);
 %put *********** Task starts at: &timenow  ************;

/* See README for more info*/

data params;
retain distribution_type avg_occurence param_1 param_2 param_3;
/*Distribution_type is the type of a continuous distrition. */
/*avg_occurence is the mean of the discrete distrition. */
/*params are the corresponding parameters in the CDF calculation in the functions. */
/*These parameters are the so called scale and shape parameters to a continuous distribution */
set params;run;

data params;
set params end=eof;
key = _N_;
If eof then call symputx('N_dist',_N_);run;


%let percentile = 0.999; /*The quantile you want*/
%let multi = 100; /*A multiplier used to extend the range of continuous distribution.*/
%let power = 20; /*The power raised on 2. Used in discrete Fourier Transform.*/

Proc iml;
*************** Foureir Transform with Poisson as Discrete Distribution***********************************;
Start Fourier(index, pdf, param, M);

	transformed_pdf = FFT(pdf);

	real_1 = exp(param[index,2] # (transformed_pdf[ ,1]-1));
	real_2 = cos(param[index,2] # transformed_pdf[ ,2]);
	compound_real = real_1 # real_2;

	imag_1 = exp(param[index,2] # (transformed_pdf[ ,1]-1));
	imag_2 = sin(param[index,2] # transformed_pdf[ ,2]);
	compound_imag = imag_1 # imag_2;

	compound_pgf = compound_real||compound_imag;
	compound_pdf = 1/M * IFFT(compound_pgf);

	sum = sum(compound_pdf);
	compound_cdf = cusum(compound_pdf);
	
	return(compound_cdf);

Finish Fourier;

************** Discretization functions and quantile calculation****************************************; 
Start High_Quantile(index, param, div_M, multi, percentile); 

**************************Lognormal Distribution***************************************;
If param[index,1]='Lognormal' then do;
	High_per = 1-(1-percentile)/param[index,2];
	High_value = Quantile('Lognormal', High_per, param[index,3], param[index,4]) # multi;
    Interval = High_value / div_M;
    y = ((0:div_M)-.5) # Interval;
    y[1,1] = 0;
    discrete_f = CDF('Lognormal', y, param[index,3], param[index,4]); 
	discrete_f = dif(discrete_f);
	discrete_f = remove(discrete_f,1);
/*	if (1 - sum(discrete_f[1:div_M])) > 1e-3 then do; print (1 - sum(discrete_f[1:div_M])) 'TOO LARGE ERROR!!'; abort; end;*/
	discrete_f[div_M] = 1 - sum(discrete_f[1:(div_M-1)]);

	compound_cdf = Fourier(index, discrete_f, param, div_M); 

 If compound_cdf[div_M-1] > percentile Then Do;

	Loss = ((1:div_M)-0.5) # Interval; 

	do i = div_M to 1 by -1;
		if compound_cdf[i] < percentile then do;
		Aggregate_Quantile = (Loss[i+1]-Loss[i])#(percentile-compound_cdf[i])/(compound_cdf[i+1]-compound_cdf[i]) + Loss[i];
		return(Aggregate_Quantile);end;
	end;
  End;

  Else Do; Print "The second CDF from the last is smaller than specified percentile. The quantile you need can't be calculated";
  		   Aggregate_Quantile = .; return(Aggregate_Quantile);
  End;

end;

*****************Burr Distribution************************************************;
If param[index,2]='Burr' then do; /*if want to limit infinite mean then add "and param_2*param_3>1"*/
	High_per = 1-(1-percentile)/param[index,2];
	High_value = param[index,3]#((1-High_per)##(-1/param[index,4])-1)##(1/param[index,5]) # multi;
    Interval = High_value / div_M;
    y = ((0:div_M)-.5) # Interval;
    y[1,1] = 0;
    discrete_f = 1 - (1/(1+(y/param[index,3])##param[index,4]))##param[index,5];
	discrete_f = dif(discrete_f);
	discrete_f = remove(discrete_f,1);
/*	if (1 - sum(discrete_f[1:div_M])) > 1e-3 then do; print 'TOO LARGE ERROR!!'; abort; end;*/
	discrete_f[div_M] = 1 - sum(discrete_f[1:(div_M-1)]);

	compound_cdf = Fourier(index, discrete_f, param, div_M);

 If compound_cdf[div_M-1] > percentile Then Do;

	Loss = ((1:div_M)-0.5) # interval; 

	do i = div_M to 1 by -1;
		if compound_cdf[i] < percentile then do;
		Aggregate_Quantile = (Loss[i+1]-Loss[i])#(percentile-compound_cdf[i])/(compound_cdf[i+1]-compound_cdf[i]) + Loss[i];
		return(Aggregate_Quantile);end;
	end;
  End;

  Else Do; Print "The second CDF from the last is smaller than specified percentile. The quantile you need can't be calculated";
  		   Aggregate_Quantile = .; return(Aggregate_Quantile);
  End;

end;

*****************Pareto Distribution************************************************;
If param[index,2]='Pareto' then do;
	High_value = param[index,3]#(((1-percentile)/param[index,2])##(-1/param[index,4])-1) # multi;
	Interval = High_value / div_M;
    y = ((0:div_M)-.5) # Interval;
    y[1,1] = 0;
	discrete_f = CDF('Pareto',y + param[index,3], param[index,4], param[index,3]);
	discrete_f = dif(discrete_f);
	discrete_f = remove(discrete_f,1);
/*	if (1 - sum(discrete_f[1:div_M])) > 1e-3 then do; print index (1 - sum(discrete_f[1:div_M])) 'TOO LARGE ERROR!!'; abort; end;*/
	discrete_f[div_M] = 1 - sum(discrete_f[1:(div_M-1)]);

	compound_cdf = Fourier(index, discrete_f, param, div_M);
	
 If compound_cdf[div_M-1] > percentile Then Do;

	Loss = ((1:div_M)-0.5) # interval; 

	do i = div_M to 1 by -1;
		if compound_cdf[i] < percentile then do;
		Aggregate_Quantile = (Loss[i+1]-Loss[i])#(percentile-compound_cdf[i])/(compound_cdf[i+1]-compound_cdf[i]) + Loss[i];
		return(Aggregate_Quantile);end;
	end;
  End;

  Else Do; Print "The second CDF from the last is smaller than specified percentile. The quantile you need can't be calculated";
  		   Aggregate_Quantile = .; return(Aggregate_Quantile);
  End;
  
End;

*****************Inverse Gaussian************************************************;
If param[index,2]='Inverse Gaussian' then do;
	High_per = 1-(1-percentile)/param[index,2];
	High_value = QUANTILE("IGAUSS", High_per, Param[index,4])#Param[index,3] # multi;
    Interval = High_value / div_M;
    y = ((0:div_M)-.5) # Interval;
    y[1,1] = 0;
    discrete_f = CDF("IGAUSS", y/param[index,3], param[index,4]);
	discrete_f = dif(discrete_f);
	discrete_f = remove(discrete_f,1);
/*	if (1 - sum(discrete_f[1:div_M])) > 1e-3 then do; print 'TOO LARGE ERROR!!'; abort; end;*/
	discrete_f[div_M] = 1 - sum(discrete_f[1:(div_M-1)]);

	compound_cdf = Fourier(index, discrete_f, param, div_M);

  If compound_cdf[div_M-1] > percentile Then Do;

	Loss = ((1:div_M)-0.5) # interval; 

	do i = div_M to 1 by -1;
		if compound_cdf[i] < percentile then do;
		Aggregate_Quantile = (Loss[i+1]-Loss[i])#(percentile-compound_cdf[i])/(compound_cdf[i+1]-compound_cdf[i]) + Loss[i];
		return(Aggregate_Quantile);end;
	end;
  End;

  Else Do; Print "The second CDF from the last is smaller than specified percentile. The quantile you need can't be calculated";
  		   Aggregate_Quantile = .; return(Aggregate_Quantile);
  End;

end;


Finish High_Quantile;


*******************************************************************************;

*****************Input**********************************************************;

percentile = &percentile;
M = 2 ## &power;
multi = &multi;
*******************************************************************************;

Use params;read all var _num_ into all_params;

Quantile_All = J(&N_dist,2,0);

DO index = 1 TO &N_dist;
	Quantile_All[index,1] = index;
	Quantile_All[index,2] = High_Quantile(index, all_params, M, multi, percentile);
END;

create Quantile_All from Quantile_All; append from Quantile_All;Print Quantile_All [format=dollar20.];
show names;

Quit;

data Quantile_All;set Quantile_All;
rename col1=primary_key
	   col2=Fourier_Quantile;
run;

proc sql;
create table Quantile_All_comb as 
select a.*,b.Fourier_Quantile format=dollar20. from
params a left join Quantile_All b
on a.key = b.primary_key;
quit;

 %let timeend=%sysfunc(time(), time.);
 %put ******* &timeend  ************;
