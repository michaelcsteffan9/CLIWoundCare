libname progress "C:\Users\haleyda\Documents\Spring 2020\BIOS 841\Work";

proc import datafile="C:\Users\haleyda\Dropbox\Consulting project\Data\analyticData.csv"
			out=progress.data
			dbms=csv
			replace;
			guessingrows=1000;
run;

%macro freq(measure=,level=,group=);
proc freq data=progress.data noprint;
tables &measure / binomial(level=&level CL=wald);
where primaryprocedure="&group";
output out=&measure binomial;
run;
%mend;

%freq(measure=ampfreesurv2yr,level=2,group=WM);
%freq(measure=majoramp2yr,level=2,group=WM);
%freq(measure=mort2yr,level=1,group=WM);
%freq(measure=mace2yr,level=2,group=WM);

data proportions;
	set ampfreesurv2yr (in=a)
		mort2yr (in=b)
		majoramp2yr (in=c)
		mace2yr (in=d);

	if a then Outcome='Two-Year Amputation Free Survival';
	if c then Outcome='Two-Year Amputation';
	if b then Outcome='Two-Year Survival';
	if d then Outcome='Two-Year MACE';

	keep outcome round _bin_ l_bin u_bin;
run;

data report;
	set proportions;
	ci='('||strip(round(l_bin,0.0001))||', '||strip(round(u_bin,0.0001))||')';
	prop=round(_bin_,0.0001);
run;

%freq(measure=ampfreesurv2yr,level=2,group=Revasc);
%freq(measure=majoramp2yr,level=2,group=Revasc);
%freq(measure=mort2yr,level=1,group=Revasc);
%freq(measure=mace2yr,level=2,group=Revasc);

data proportions;
	set ampfreesurv2yr (in=a)
		mort2yr (in=b)
		majoramp2yr (in=c)
		mace2yr (in=d);

	if a then Outcome='Two-Year Amputation Free Survival';
	if c then Outcome='Two-Year Amputation';
	if b then Outcome='Two-Year Survival';
	if d then Outcome='Two-Year MACE';

	keep outcome round _bin_ l_bin u_bin;
run;

data report_2;
	set proportions;
	ci='('||strip(round(l_bin,0.0001))||', '||strip(round(u_bin,0.0001))||')';
	prop=round(_bin_,0.0001);
run;


ods rtf file="C:\Users\haleyda\Documents\Spring 2020\BIOS 841\Work\aim1.rtf";

proc sgplot data=proportions;
	label _bin_ = 'Proportion';
	scatter y=outcome x=_bin_ / xerrorlower=l_bin xerrorupper=u_bin markerattrs=(symbol=circlefilled) errorbarattrs=(color=viyg);
	yaxis discreteorder=data offsetmax=0.1 display=(nolabel);
	xaxis grid values=(0 to 1 by 0.1) valueshint;
run;

title 'Wound Management Group';
proc report data=report nowd
	style(column header)=[background=white];
	columns outcome prop ci;
	define outcome / display 'Outcome';
	define prop / display 'Proportion';
	define ci / display '95% Wald Confidence Interval';
run;

title 'Revascularization Group';
proc report data=report_2 nowd
	style(column header)=[background=white];
	columns outcome prop ci;
	define outcome / display 'Outcome';
	define prop / display 'Proportion';
	define ci / display '95% Wald Confidence Interval';
run;
title;

ods rtf close;
