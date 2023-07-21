/*use with imagej 1.53 or up (done with imagej 1.53f)
 *This Macro finds the center coordinates of a given cell in a crossbow-shaped pattern of images 512x512 pixels. 
 *Additionally it selects a focused Z and rotates the outer edge of the cell to match the center
 *To find the center of the cell, the cell must be positioned in the 4th quartile
 *After the initial rotation to the 4th quartile a new rotation of 90 degrees will be initated 
 *The center of the cell's edge is determined by 3 rectangles which measure the total coverage 
 *The angle will be selected which gave the max coverage of the three rectangles.
 *
 *
 *
 * Scriped by Daniel SAMPAIO GONCALVES (daniels_94@hotmail.fr) 01/12/2020 with ImageJ version 1.53
 ***********************************************************************************************************/
requires("1.53"); 
setBatchMode(true);
print("\\Clear");
start = getTime();

setOption("ExpandableArrays", true); //is needed to expand arrays which will be used down below


Area1=newArray;
Area2=newArray;
Area3=newArray;


dirIn=getDirectory("Folder with actin Z stack");


Probability=getDir("Probabilitymaps");

filesTmp  = getFileList(dirIn);
files     = filterFileNames(filesTmp,".tif",true);


print("\\Clear");


  for(i=0; i<files.length;i++)
{     
	  
	open(dirIn+files[i]);


  //Selects focused slice
 	print("N"+i+1+"= "+files[i]);
  	run("Find focused slices", "select=80 variance=0.000");
	NewImageTitle=getTitle();
	selectWindow(NewImageTitle);
	close();
	z0=getSliceNumber();

	setSlice(z0);
	//run("Duplicate...", "title=focus.tif"); //for debug purposes below, focus.tif will be created by the MakeMask() function
	
	run("Duplicate...", "title=focus1.tif");
	MakeMask(Probability,files,i);

   
   //selectWindow("focus.tif"); //for debug purposes
   //setAutoThreshold();
   //run("Convert to Mask");
   //run("Fill Holes");
  
   
//calculates the center coordinates
 
selectWindow("focus1.tif");
run("Set Scale...", "distance=0 known=0 unit=pixel");
setAutoThreshold();
run("Find Edges");
run("Set Measurements...", "center redirect=None decimal=10");

   print("here z0 = ", z0);
   print( "nResults = ", nResults );
   print( "i = ", i );

   setResult("X", i, getValue("XM"));
   setResult("Y", i, getValue("YM"));
   setResult("Z", i, z0);

	
selectWindow("focus.tif");

//translates the image into the center coordinates

run("Fill Holes");
getDimensions(width, height, channels, slices, frames);
x0=parseInt(getResult("X",i));
y0=parseInt(getResult("Y",i));
print("y0="+y0);
print("x0="+x0);
run("Translate...", "x="+(width/2-x0)+" y="+(height/2-y0)+" interpolation=None");


selectWindow("focus.tif");

//this portion will initiate the initial rotation angle based on images sets
Current=ReplaceToNumber(files,i);
		
		if (i-1<0) { 
		Previous=Current;
		AA=InitialRotationCheck();
		}
		
		else {Previous=ReplaceToNumber(files,i-1);}
		
		
		if (Current!=Previous) {AA=InitialRotationCheck();}

		
		BB=AA;
		A=abs(AA)*PI/180;

run("Rotate... ", "angle="+BB+" grid=1 interpolation=Bilinear"); //rotates to the initial selected angle
print("First Rotation using a ZAngle of: " +(A)+" Rad "+" or "+BB+" Degree");

max3=0;
y1=390;
ti=-1;
while (max3<70) {

y1=y1-10;	
ti=ti+1;
run("Duplicate...", "title=rotate.tif");

for (Rot = 0; Rot < 90; Rot++) //loop that rotates 1 degree each loop cycle for 90 degrees and adds all the coverage values of the three rectangles to 3 arrays
{  
	makeRectangle(216, y1, 80, 8);
	MeanCenter1=getValue("%Area");
	roiManager("add");
	roiManager("Show All");
	makeRectangle(216, y1+30, 80, 8);
	MeanCenter2=getValue("%Area");
	makeRectangle(216, y1+60, 80, 8);
	MeanCenter3=getValue("%Area");
	run("Select None");
	

	
	 Area1=addToArray(parseInt(toString(MeanCenter1,0)),Area1,  Rot);
	 Area2=addToArray(parseInt(toString(MeanCenter2,0)),Area1,  Rot);
	 Area3=addToArray(parseInt(toString(MeanCenter3,0)),Area1,  Rot);

	 
	 run("Rotate... ", "angle="+1+" grid=1 interpolation=Bicubic"); 
}
roiManager("Delete");
close();

//get stats of the array to print the maximum of the array
Array.getStatistics(Area1, min, max1, mean, stdDev);
Array.getStatistics(Area2, min, max2, mean, stdDev);
Array.getStatistics(Area3, min, max3, mean, stdDev);
print("Maximum % cell coverage at Positions : A1: "+max1+"% A2: "+max2+"% A3: "+max3+"%; moved up= "+ti+" times");

//if (max3<60) {waitForUser;}  //for debugging purposes

}

pos1=Array.findMaxima(Area1, 1); //returns the position of the maximum found coverage)
pos2=Array.findMaxima(Area2, 1);
pos3=Array.findMaxima(Area3, 1);


print("Angle to correct to reach maximum coverage Pos1: "+pos1[0]+" Pos2; " +pos2[0]+" Pos3: "+pos3[0]);

BB=(pos1[0]+pos2[0]+pos3[0])/3; //averages the position of the 3 rectangles , or here the angle needed to get max coverage





A=abs(AA+BB)*PI/180;

	
print("Second Rotation using a ZAngle of: " +(BB)*PI/180+" Rad "+" or "+ abs(BB)+" Degree");




setResult("ZAngle",i, 2*PI+(-A+PI)); //sets the final angle into the result table 

print("New ZAngle: " + (A)+ " Rad or "+ (180+BB)+" Degrees"+"\nAdditionally 2*PI+(-[ZAngle] +PI)="+2*PI+(-A+PI)+" was calculated to correct for the 2nd Macro\n");



close();


updateResults(); //updates result table


//to debug the final output open this line
/* 
selectWindow("TEST.tif");
print((A*180)/PI);
waitForUser("TESTING");
run("Rotate... ", "angle="+(A*180)/PI+" grid=1 interpolation=Bicubic");
run("Grid...", "grid=Lines area=8000 color=Cyan center");
waitForUser("ROTATION");
*/

while(nImages>0){
   close();
}
run("Collect Garbage"); //empties imagej memory cash, useful especially for MacOs when opening a lot of files
}


saveAs("Results", dirIn+"../Results.csv"); //save results in the folder before the image folder



//function

function filterFileNames(filelist,extension,ignorecase)
 //filters files, only opens tifs
{
	  tmpfiles=newArray(filelist.length);
	  count=0;
	  pattern=extension;
	  
	  for(i=0;i<filelist.length;i++)
	  {
		    strin=filelist[i];
		    if (ignorecase)
		    {
			      strin=toLowerCase(strin);
			      pattern=toLowerCase(pattern);
		    }
		    
		    if (endsWith(strin,pattern))
		    {
			      tmpfiles[count]=filelist[i];
			      count++;
		    }
	  }
	  
	  filteredList=newArray(count);
	  for (i=0;i<count;i++)	filteredList[i]=tmpfiles[i];
	  
	  return filteredList;
}



function addToArray(value, array, position) //function to add values to a given array
{
    if (position<lengthOf(array)) {
        array[position]=value;
    } else {
        temparray=newArray(position+1);
        for (i=0; i<lengthOf(array); i++) {
            temparray[i]=array[i];
        }
        temparray[position]=value;
        array=temparray;
    }
    return array;
}

function InitialRotationCheck() //function to open, check and set the initial rotation angle
{
run("Duplicate...", "title=TestRotation.tif");
setBatchMode("show");
Dialog.create("Chose Angle");
Dialog.addMessage("The edge of the cell needs\nto be on the 4th Quartile\nAvoid touching the ROI");
Dialog.addNumber("Angle:", 180);
Dialog.show();

ANGLE=Dialog.getNumber();

run("Rotate... ", "angle="+ANGLE+" grid=1 interpolation=Bilinear");
run("Grid...", "grid=Lines area=8000 color=Cyan center");
makeRectangle(216, 390, 80, 8);
roiManager("add");
makeRectangle(216, 420, 80, 8);
roiManager("add");
makeRectangle(216, 450, 80, 8);
roiManager("add");
run("Select None");
roiManager("show all");

ITEMS=newArray("Yes","No");

Dialog.create("Is the Angle correct?");
Dialog.addMessage("The edge of the cell needs\nto be on the 4th Quartile\nAvoid touching the ROI");
Dialog.addRadioButtonGroup("Is the Angle Correct?", ITEMS, 1, 2, "Yes");
Dialog.show();
Response=Dialog.getRadioButton();
close();
roiManager("delete");
while (Response=="No") { 
	
run("Duplicate...", "title=TestRotation.tif");
setBatchMode("show");
Dialog.create("Chose Angle");
Dialog.addMessage("The edge of the cell needs\nto be on the 4th Quartile\nAvoid touching the ROI");
Dialog.addNumber("Angle:", 0);
Dialog.show();

ANGLE=Dialog.getNumber();

run("Rotate... ", "angle="+ANGLE+" grid=1 interpolation=Bilinear");
run("Grid...", "grid=Lines area=8000 color=Cyan center");

makeRectangle(219, 370, 42, 4);
roiManager("add");
makeRectangle(219, 375, 42, 4);
roiManager("add");
makeRectangle(219, 380, 42, 4);
roiManager("add");
run("Select None");
roiManager("show all");

Dialog.create("Is the Angle correct?");
Dialog.addMessage("The edge of the cell needs\nto be on the 4th Quartile\nAvoid touching the ROI");
Dialog.addRadioButtonGroup("Is the Angle Correct?", ITEMS, 1, 2, "Yes");
Dialog.show();
Response=Dialog.getRadioButton();
close();
roiManager("delete");
}
setBatchMode("hide");
return ANGLE ;
}

function ReplaceToNumber(files,i)//Function to extract the number of the different image sets; names should be in this scheme E20-42_WM1862_I_01_01_R3D_D3D_W525.tif. Numbers and letters are case sensitive!
{
	

		a=replace(files[i], "[A-Z]+_[0-9][0-9]_[0-9][0-9]_[A-Z][0-9][A-Z]_[A-Z][0-9][A-Z].*", "");
		b=replace(files[i], "[A-Z][0-9][0-9]-[0-9][0-9].*_[A-Z]+_[0-9][0-9]_[0-9][0-9]_", "");
		c=replace(replace(files[i], b, ""), a, "");
		x=replace(c, "_[0-9][0-9]_[0-9][0-9]_", "");
		//c=replace(replace(c, x+"$", ""), "_", "");

		return x;
}
function MakeMask(path,files,i) { 
	

open(path+files[i]);
/*
setAutoThreshold("Otsu dark"); 
*/
setThreshold(0.5, 1);

run("Convert to Mask", "method=Otsu background=Dark");
run("Duplicate...", "  channels=1 title=focus.tif");	
selectWindow("focus.tif");
run("Fill Holes");

}


print("**Processing finished, took "+((getTime()-start)/1000)+" second(s) **");