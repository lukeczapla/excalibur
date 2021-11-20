#ifndef _VIS_H
#define _VIS_H

#include "gl.h"
#include "../util/socket.h"
#include "../util/thread.h"

#define MAX_TYPES 100


void ViewOrtho(int width, int height, float factor);

extern int StringToInt(char * str, int base, char *file, int lineno);
extern double StringToDouble(char *str, char *file, int lineno);
extern int IsDelim(char test, char *delimiters, int ndelim);
extern int TokenizeString(char *str, char *delimiters, char **pointers, int maxptrs);
extern float ran1(long *idum);


/*
	*****************************************************************************************************************
	Simple GL window class derived from GLWindow. Has complete user interaction via mouse and keyboard, and an Apple
	specific vertical blank sync.
	*****************************************************************************************************************
*/
class vis : public GLWindow {

	public:
		GLfloat lastx, lasty;
		GLint display_timer_id, move_timer_id, info_timer_id, framecount, target_fps;

		// scene info and general draw variables
		GLint scene_dlist;
		bool data_invalidated, view_invalidated, vblank_sync;
				
		// cell info
		double cell[3];
		GLfloat cell_col[3];
		
		// site info
		int nsites, *types;
		bool show_types[MAX_TYPES];
		double *r, *rad;
		
		long ran1_seed;

		// bead info to draw sites
//		GLint bead_dlist;
		GLfloat bead_amb[4*MAX_TYPES];
		GLfloat bead_diff[4*MAX_TYPES];
		GLfloat bead_spec[4*MAX_TYPES];
		GLfloat bead_emit[4*MAX_TYPES];
		GLfloat bead_shiny[MAX_TYPES];
		
		char frame_string[50];
		
		vis(int _target_fps, int argc, char **argv);
		~vis();
		void Timer(GLint value);
		void Display();
		void MouseButton(GLint button, GLint state, GLint x, GLint y);
		void MouseMotion(GLint x, GLint y);
		void KeyDown(GLint key, GLint x, GLint y);
		void Reshape(GLint width, GLint height);
		
		void setLightPos(GLfloat x, GLfloat y, GLfloat z);	
		void Update(int _n, double * _r , int * _t, double * _cell, double *_rad);
};


#endif

