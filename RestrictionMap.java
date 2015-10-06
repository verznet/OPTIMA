/*  Type: RestrictionMap

    Set score limits and restriction map information.

    
    --
    OPTIMA v.f-1.3 -- 6 October 2015
    Copyright (C) Davide Verzotto and Niranjan Nagarajan

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License 2.1 as published by the Free Software Foundation.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License 2.1 along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301, USA.

 */
 
 
package OPTIMA;

import java.text.DecimalFormat;

public class RestrictionMap
{  
   public static int MINIMUMSCORE = -10000000;
   public static int maxNoOfMatches = 30000; //500;
   public static int ENDOFCHROMOSOMES = 11000000;

   int noFragments;
   Fragment[] fragments;
   int index = 0;
   String id = "";
   int sizeMapInBases = 0;
   

   public RestrictionMap(Fragment[] setOfFragments, int index, String id)
   { fragments = setOfFragments;
     noFragments = fragments.length;
   }
   
   public RestrictionMap()
   { noFragments = 0;
     fragments = new Fragment[30];
   }
   
   public Fragment get(int i)
   { return fragments[i];
   }
   
   public int size()
   { return noFragments;
   }
   
   static final DecimalFormat numberFormat = new DecimalFormat("#.000");
	
   public String toString() {
		String str = "";
		str = this.id + "\n" + "\tKpnI" + "\tK"; // + "\tnoOfFragments:" + this.noFragments + "\n";
		for (Fragment frag : this.fragments)
		{ str = str + "\t" + numberFormat.format((frag.size / 1000)); }
		
		str = str + "\n"; 
		return str;
   }

}
