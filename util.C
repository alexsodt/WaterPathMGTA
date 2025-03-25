// by alex sodt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#include <stdarg.h>
int readNDoubles( char *buffer, double *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%lf", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}



void getLine( FILE *theFile, char *theBuffer )
{
        int i = 0;

        while( !feof(theFile) )
        {
                char tc = fgetc(theFile);

                if( tc != '\n' && i < 1000000 )
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n' && i >= 1000000 )
		{
		}
		else
                        break;
        }

        theBuffer[i] = '\0';
}

int goToField( const char * buffer, int f )
{
	int ret = -1;
	int cf = 0;
	
	int cp = 0;
	while( buffer[cp] && (buffer[cp] == ' ' || buffer[cp] == '\t') )
		cp++;

	if( !buffer[cp] ) return -1;


	while( cf < f )
	{
		while( buffer[cp] && (buffer[cp] != ' ' && buffer[cp] != '\t') )
			cp++;
		
		if( !buffer[cp] ) return -1;

		while( buffer[cp] && (buffer[cp] == ' ' || buffer[cp] == '\t') )
			cp++;

		if( !buffer[cp] ) return -1;

		cf++;
	}
		

	return cp;
}

int readNInts( char *buffer, int *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
        int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {   
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%d", vals+x );

                if( nr != 1 ) 
                        return nread;

                nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }   

        return nread;
}

int my_isnan( double a )
{
	if( !(a > 0 || a < 1) )
		return 1;
	return 0;
}


static FILE *default_fp = 0;

void set_default_stream(FILE *fp)
{
    default_fp = fp;
}

int mprintf(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    if (default_fp == 0)
        default_fp = stdout;

    int rv = vfprintf(default_fp, fmt, args);

    va_end(args);
    return(rv);
 }

int decodeString( char *buf, char **out, char **res_out, int nmax )
{
	int done = 0;

	int ns = 0;
	int t = 0;
	int lt = 0;
	int lim = 256;
	char cur[lim+1];
	char res_cur[lim+1];

	while( *(buf+t) )
	{
		while( *(buf+t) == '+' ) t++;

		int got_res = 0;

		while( *(buf+t) && *(buf+t) != '+' && *(buf+t) != '_' )
		{
			if( lt < lim )
			{
				cur[lt] = *(buf+t);	
				lt++;
			}

			t++;
		}
		
		cur[lt] = '\0';

		int res_lt = lt;

		if( *(buf+t) == '_' )
 		{
			t += 1;
			got_res = 1;

			strcpy( res_cur, cur );

			lt = 0;

			while( *(buf+t) && *(buf+t) != '+' )
			{
				if( lt < lim )
				{
					cur[lt] = *(buf+t);	
					lt++;
				}
	
				t++;
			}
		
			cur[lt] = '\0';
		}

 
		if( lt > 0 && ns < nmax )
		{
			out[ns] = (char *)malloc( sizeof(char) * (1 + lt ) );
			strcpy( out[ns], cur );

			if( got_res ) 
			{						
				res_out[ns] = (char *)malloc( sizeof(char) * (1 + res_lt ) );
				strcpy( res_out[ns], res_cur );
			}
			else
			{
				res_out[ns] = (char *)malloc( sizeof(char) * (1 + strlen("any") ) );
				strcpy( res_out[ns], "any" );
			}
			ns++;
		}
			
		lt = 0;

	}

	return ns;
}
/*
int decodeString( char *buf, char **out, int nmax )
{
	int done = 0;

	int ns = 0;
	int t = 0;
	int lt = 0;
	int lim = 256;
	char cur[lim+1];

	while( *(buf+t) )
	{
		while( *(buf+t) == '+' ) t++;

		while( *(buf+t) && *(buf+t) != '+' )
		{
			if( lt < lim )
			{
				cur[lt] = *(buf+t);	
				lt++;
			}

			t++;
		} 		

		cur[lt] = '\0';
 
		if( lt > 0 && ns < nmax )
		{
			out[ns] = (char *)malloc( sizeof(char) * (1 + lt ) );
			strcpy( out[ns], cur );
			ns++;
		}
			
		lt = 0;

	}

	return ns;
}
*/
void print5( int val, char *str )
{
        if( val < 10 )
                sprintf(str, "0000%d", val );
        else if( val < 100 )
                sprintf(str, "000%d", val );
        else if( val < 1000 )
                sprintf(str, "00%d", val );
        else if( val < 10000 )
                sprintf(str, "0%d", val );
        else
                sprintf(str, "%d", val );
}

int read_sequence( char **seq, FILE *theFile )
{
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	
	getLine( theFile, buffer );

	if( buffer[0] != '>' )
		return 0;

	int fpos = ftell(theFile);
	getLine( theFile, buffer );

	do { 
		if( ! (*seq) )
		{
			*seq = (char *)malloc( sizeof(char) * (1 + strlen(buffer) ) );
			strcpy( *seq, buffer );
		}
		else
		{
			int cur_len = strlen(*seq);
			*seq = (char *)realloc( *seq, sizeof(char) * ( 1 + cur_len + strlen(buffer) ) );
			strcpy( *seq+cur_len, buffer );
		}		

		fpos = ftell(theFile);
		getLine( theFile, buffer );
	} while( !feof(theFile) && buffer[0] != '>' );

	if( buffer[0] == '>' )
		fseek( theFile, fpos, SEEK_SET );

	free(buffer);

	return strlen(*seq); 
}
