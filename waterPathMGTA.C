#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "alignSet.h"
#include "mutil.h"

#define WATER_BULK	0
#define WATER_UPPER	1
#define WATER_TM	2
#define WATER_LOWER	3

static double water_upper_boundary;
static double water_lower_boundary;
static double cutoff=3.4;
static double cutoff_SC=5.0;
int main( int argc_in, char **argv_in )
{
	int argc = 0;
	char *argv[argc_in];
	int nres_seek = 0;
	int max_res_seek = 10;
	int res_seek[max_res_seek];
	double boundary = 20;
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	int doB=0;
	int minw = 2;
	int cterm = 0;
	int start = -1;
	int stop = -1;
	int totalPath = 0;	


	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2 ) )
		{
			if( !strncasecmp( argv_in[c], "--boundary=", strlen("--boundary=") ) ) 
				boundary = atof(argv_in[c]+strlen("--boundary="));
			else if( !strncasecmp( argv_in[c], "--cutoff=", strlen("--cutoff=") ) ) 
				cutoff = atof(argv_in[c]+strlen("--cutoff="));
			else if( !strncasecmp( argv_in[c], "--minw=", strlen("--minw=") ) ) 
				minw = atoi(argv_in[c]+strlen("--minw="));
			else if( !strncasecmp( argv_in[c], "--start=", strlen("--start=") ) ) 
				start = atoi(argv_in[c]+strlen("--start="));
			else if( !strncasecmp( argv_in[c], "--stop=", strlen("--stop=") ) ) 
				stop = atoi(argv_in[c]+strlen("--stop="));
			else if( !strncasecmp( argv_in[c], "--seek=", 7 )  )
			{
				if( nres_seek < max_res_seek )
				{
					res_seek[nres_seek] = atoi(argv_in[c]+7);
					nres_seek++;
				}
			}
			else if( !strcasecmp( argv_in[c], "--B" )  )
				doB = 1;
			else if( !strcasecmp( argv_in[c], "--cterm" )  )
				cterm = 1;
			else if( !strcasecmp( argv_in[c], "--totalPath" )  )
				totalPath = 1;
			else 
			{
				printf("Couldn't interpret flag '%s'.\n", argv_in[c] );
				exit(1);
			}
		}
		else
		{
			argv[argc] = argv_in[c];
			argc++;
		}
	}
	if( argc < 3 )
	{
		printf("Syntax: waterPathMGTA psf dcd \n");	
		return 0;
	}
	water_upper_boundary = boundary;
	water_lower_boundary = -boundary;

	FILE *psfFile = fopen(argv[1],"r");
	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );
		
	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	int nat = curNAtoms();
	int init_done = 0;

	int nwat = 0;
	int nwatSpace = 10;
	int *wat = (int*)malloc( sizeof(int) * nat );
	int mg_atom =-1;
	int res_atom[max_res_seek];

	double nupper = 0;
	double nlower = 0;
	double ntotal = 0;

	double *meets = (double *)malloc( sizeof(double) * max_res_seek * max_res_seek );
	memset( meets, 0, sizeof(double) * max_res_seek * max_res_seek );

	int ftot = 0;
	
	// for finding contacts.
	//
	struct binnit
	{
		double *pos;
		int *id;
		int *link;
		int npos;
		int nposSpace;
	};
	binnit *bins = NULL;
	int nx,ny,nz;
	double tbinw = 7.0;

	for( int c = 2; c < argc; c++ )
	{
		FILE *dcdFile = fopen( argv[c], "r");
	
		if( !dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[3] );
			exit(1);
		}
		readDCDHeader(dcdFile);
	
		setAligned();
	
		int nframes = curNFrames();
	
		for( int f = 0; f < nframes; f++, ftot++)
		{
			if( stop >= 0 && ftot > stop )
				break;

			loadFrame( dcdFile, at );
		
			double last_PBC[3];
			double cur_PBC[3];
			double La,Lb,Lc,alpha,beta,gamma;	
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
		
			if( !DCDsuccess() )
			{
				printf("Could not load intended frames\n");
				exit(1);
			}
	
			if( start >= 0 && ftot < start )
			{
				for( int a = 0; a < nat; a++ )
					at[a].zap();
				continue;
			}
		
			if( ! init_done )
			{
				nx = La/tbinw;
				ny = Lb/tbinw;
				nz = Lc/tbinw;

				bins = (binnit *)malloc( sizeof(binnit) * nx * ny * nz );

				for( int x = 0; x < nx*ny*nz; x++ )
				{
					bins[x].nposSpace = 20;
					bins[x].id = (int*)malloc( sizeof(int) * bins[x].nposSpace );
					bins[x].npos = 0; 
				}

				for( int a =0;a<nat;a++)
				{
					if(!strcasecmp(at[a].resname,"TIP3")&&at[a].atname[0]=='O')
					{
						if( nwatSpace == nwat )
						{
							nwatSpace *= 2;
							wat = (int *)realloc(wat,sizeof(int)*nwatSpace);
						}
						wat[nwat++]=a;
					}
					if( !strcasecmp(at[a].atname,"MG") && ((!doB && !strcasecmp(at[a].segid, "HETA" )) || (doB && !strcasecmp(at[a].segid,"HETB"))) )
						mg_atom = a;
					if( !strcasecmp(at[a].atname,"MG") && ((!doB && !strcasecmp(at[a].segid, "MGA" ) && at[a].res == 1) || (doB && !strcasecmp(at[a].segid,"MGB") && at[a].res == 1 )) )
						mg_atom = a;

					for( int r = 0; r < nres_seek; r++ )
					{
						if( at[a].res == res_seek[r] && ((!doB && !strcasecmp(at[a].segid,"PROA")) || (doB && !strcasecmp(at[a].segid,"PROB"))) )
						{
							if( !strcasecmp( at[a].resname, "GLU" ) && !strcasecmp( at[a].atname, "CD") )
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "ASP" ) && !strcasecmp( at[a].atname, "CG") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "SER" ) && !strcasecmp( at[a].atname, "OG") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "SER" ) && !strcasecmp( at[a].atname, "OG") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "THR" ) && !strcasecmp( at[a].atname, "OG1") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "TYR" ) && !strcasecmp( at[a].atname, "OH") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "ARG" ) && !strcasecmp( at[a].atname, "CZ") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "LYS" ) && !strcasecmp( at[a].atname, "NZ") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "ASN" ) && !strcasecmp( at[a].atname, "CG") )	
								res_atom[r] = a;
							else if(  !strcasecmp( at[a].resname, "GLN" ) && !strcasecmp( at[a].atname, "CD") )	
								res_atom[r] = a;
						} 
					}
	
				}
				if( mg_atom == -1 )
				{
					printf("Couldn't find the right Mg ion to use.\n");
					exit(1);
				}
				init_done=1;
			}
			
			for( int b = 0; b < nx*ny*nz; b++ )
				bins[b].npos = 0;
			
			for( int aw = 0; aw < nwat; aw++ )
			{
				int a = wat[aw];
				
				double x = at[a].x/La;
				double y = at[a].y/Lb;
				double z = at[a].z/Lc;
				while( x < 0 ) x += 1;
				while( x >= 1.0 ) x -= 1;
				while( y < 0 ) y += 1;
				while( y >= 1.0 ) y -= 1;
				while( z < 0 ) z += 1;
				while( z >= 1.0 ) z -= 1;

				int bx = nx * (x);
				int by = ny * (y);
				int bz = nz * (z);
				if( bx < 0 ) bx = 0; if ( bx >= nx ) bx = nx-1;
				if( by < 0 ) by = 0; if ( by >= ny ) by = ny-1;
				if( bz < 0 ) bz = 0; if ( bz >= nz ) bz = nz-1;

				int tbin = bz+nz*(by+ny*bx);
		
				if( bins[tbin].nposSpace == bins[tbin].npos )
				{
					bins[tbin].nposSpace *= 2;
					bins[tbin].id = (int *)realloc( bins[tbin].id, sizeof(int) *  bins[tbin].nposSpace );
				}
		
				bins[tbin].id[bins[tbin].npos] = aw;		
				bins[tbin].npos += 1;
			}
		
			int *wat_all    = (int*)malloc( sizeof(int) * nwat );
			int *wat_region = (int*)malloc(sizeof(int)*nwat);
			int nwatr=0;
		
			for( int aw = 0; aw < nwat; aw++ )
			{
				int a = wat[aw];
					
				double z = at[a].z;

				if( z < water_lower_boundary - 10 ) 
					wat_all[aw] = WATER_BULK;
				else if( z < water_lower_boundary )
					wat_all[aw] = WATER_LOWER;
				else if( z < water_upper_boundary )
					wat_all[aw] = WATER_TM;
				else if( z < water_upper_boundary + 10 )
					wat_all[aw] = WATER_UPPER;
				else
					wat_all[aw] = WATER_BULK;
					
				if( at[a].z > water_lower_boundary && at[a].z < water_upper_boundary )
					wat_region[nwatr++]=a;
			}

			int done = 0;
			int max_borders = 10;

			// utility mark, for example, for C-terminus			
			int *mark = (int *)malloc( sizeof(int) * nwatr );
			memset( mark, 0, sizeof(int) * nwatr );

			int *borders = (int *)malloc( sizeof(int) * max_borders * nwatr );
			int *nborders = (int *)malloc( sizeof(int) * nwatr );
			memset( nborders, 0, sizeof(int) * nwatr );
			int *depth[1+max_res_seek];

			for( int i = 0; i < 1 + nres_seek; i++ )
				depth[i] = (int*)malloc(sizeof(int)*nwatr);


			int unset = 10000;

			for( int aw = 0; aw < nwatr; aw++ )
			{
				int a = wat_region[aw];

				for( int bw = aw+1; bw < nwatr; bw++ )
				{
					int a2 = wat_region[bw];
					double dr[3] = { at[a].x-at[a2].x,at[a].y-at[a2].y,at[a].z-at[a2].z};
					double r = normalize(dr);

					if( r < cutoff )
					{
						borders[aw*max_borders+nborders[aw]] = bw;
						borders[bw*max_borders+nborders[bw]] = aw;
						nborders[aw]++;
						nborders[bw]++;
					}							
				}					
			}
				
			for( int rx = 0; rx < 1 + nres_seek; rx++ )
			{
				for( int aw = 0; aw < nwatr; aw++ )
					depth[rx][aw] = unset;

				if( rx == 0 ) // the Magnesium ion
				{
					for( int aw = 0; aw < nwatr; aw++ )
					{
						int a = wat_region[aw];
						double dr[3] = { at[a].x-at[mg_atom].x,at[a].y-at[mg_atom].y,at[a].z-at[mg_atom].z};
						double r = normalize(dr);
		
						if( r < cutoff )
							depth[rx][aw] = 0;	
					}
				}
				else
				{
					for( int aw = 0; aw < nwatr; aw++ )
					{
						int a = wat_region[aw];
						double dr[3] = { at[a].x-at[res_atom[rx-1]].x,at[a].y-at[res_atom[rx-1]].y,at[a].z-at[res_atom[rx-1]].z};
						double r = normalize(dr);
		
						if( r < cutoff_SC )
							depth[rx][aw] = 0;	
					}
				}
				done = 0;
	
				while( !done )
				{
					done = 1;
					for( int aw = 0; aw < nwatr; aw++ )
					{
						if( depth[rx][aw] == unset  )
							continue;
	
						if( nborders[aw] < minw )
							continue;

						for( int bx = 0; bx < nborders[aw]; bx++ )
						{
							if( nborders[borders[aw*max_borders+bx]] < minw )
								continue;

							if(  depth[rx][borders[aw*max_borders+bx]] > depth[rx][aw]+1 )
							{
								depth[rx][borders[aw*max_borders+bx]] = depth[rx][aw]+1;
								done = 0;
							} 
						}
					}
				}
			}
			
			double l_meets[(1+nres_seek)*(1+nres_seek)];

			for( int r1 = 0; r1 <= nres_seek; r1++ ) 
			for( int r2 = 0; r2 <= nres_seek; r2++ ) 
			{
				l_meets[r1*(1+nres_seek)+r2] = 0;

				for( int aw = 0; aw < nwatr; aw++ )
				{
					if( depth[r1][aw] != unset && depth[r2][aw] != unset )
					{
						l_meets[r1*(1+nres_seek)+r2] = 1;
						break;
					}	
				}
			}
	
			for( int r1 = 0; r1 <= nres_seek; r1++ ) 
			for( int r2 = 0; r2 <= nres_seek; r2++ ) 
				meets[r1*(1+nres_seek)+r2] += l_meets[r1*(1+nres_seek)+r2];

			int max_depth = 20;
			int path_depth[max_depth];
			memset( path_depth, 0, sizeof(int) * max_depth );
		
			int meets_with_upper_surface = 0;
			int meets_with_lower_surface = 0;
	
			// is there a water path to +/- BORDER?
			for( int aw = 0; aw < nwatr; aw++ )
			{
				if( depth[0][aw] >= 0 && depth[0][aw] < max_depth ) 
				{
					int a = wat_region[aw];

					if( at[a].z > water_upper_boundary - cutoff )
					{
						meets_with_upper_surface = 1;
					}
					if( at[a].z < water_lower_boundary + cutoff )
					{
						meets_with_lower_surface = 1;
					}
				}
                                if( depth[0][aw] >= 0 && depth[0][aw] < max_depth ) 
                                        path_depth[depth[0][aw]] += 1;
			}

			printf("frame %d.\n", f);
			for( int b = 0; b < max_depth; b++ )
			{
				printf("%d %d\n", b, path_depth[b] );
			}
			free(wat_region);

			if( meets_with_upper_surface ) { nupper += 1; printf("UPPER %s frame %d\n", argv[c], f); }
			if( meets_with_lower_surface ) { nlower += 1; printf("LOWER %s frame %d\n", argv[c], f); }
			ntotal++;

			for( int a = 0; a < nat; a++ )
				at[a].zap();
		}
		fclose(dcdFile);
	}

	printf("Fraction meeting upper %lf lower %lf\n", nupper / ntotal, nlower / ntotal );

	if( nres_seek > 0 )
	{
		printf("Res seek pairing:\n");
		for( int r1 = 0; r1 <= nres_seek; r1++ )
		{
			if( r1 == 0 )
				printf("Mg " );
			else
				printf("%d ", at[res_atom[r1-1]].res );
			for( int r2 = 0; r2 <= nres_seek; r2++ )
				printf("%lf ", meets[r1*(1+nres_seek)+r2] / ntotal );
			printf("\n");
		}

	}
}




