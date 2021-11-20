#ifndef _VIS_GL_H
#define _VIS_GL_H

#define ORTHOGRAPHIC 0

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#ifdef __APPLE__
	#include <OpenGL/OpenGL.h>
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
	#include <GLUT/glut.h>
#else
	#include <GL/gl.h>
	#include <GL/glu.h>
	#include <GL/glut.h>
#endif

/*
	OVERVIEW:
	
	Standard _glut[X]() c functions used as hooks into the GLUT callbacks;
	these functions call the _glut[X]() methods of _currentGLWindow, which in turn call
	the (overwritten) [X]() methods of that instance. Subclass GLWindow and overwrite the virtual methods
	for the desired behaviour.
	
	_currentGLWindow is a global pointer to the active GLWindow. You can change this using SetCurrentGLWindow(),
	but any previous GLWindow will lose the GLUT messages.
	
	Worth knowing:
	
	- Creation of a GLWindow instance automatically sets it to _currentGLWindow.
	- The "keys" arrays in the GLWindow class stores the current state of the keyboard so you can check
		this array where desired, instead of relying on responding to keyboard events. This array is altered
		automatically in the _glut[X]() methods, so subclased GLWindows maintain this behaviour with no user
		code required.
	- The "windowdims", "mousebuttons" and "mousepos" arrays behave in a similar manner to the above "keys" array.
	- Although it's possible to overwrite the _glut[X]() methods in GLWindow, this is not recommended, as you can lose the
	    automatic updates to the arrays mentioned above, and/or the timers depending on what you change.
	    Try to stick to overwriting the [X]() methods!
	
	To do:
	
	- Change _currentGLWindow to an array of pointers for multiple updaing instances of GLWindow.
	- Fix timers; allow the switching off of timers ( set the appropriate entry in the "timers" array )
	- VBL sync with sentinels ( varies from OS to OS ).
*/

class v3;
class camera;
class light;
class GLWindow;
extern GLWindow * _currentGLWindow;

extern void ErrorMessage(const char *calledfrom, const int linenumber, const char *msg);
extern void SetCurrentGLWindow(GLWindow *glw);
extern void _glutTimerHook(GLint value);
extern void _glutReshapeHook(GLint width, GLint height);
extern void _glutDisplayHook();
extern void _glutMouseButtonHook(GLint button, GLint state, GLint x, GLint y);
extern void _glutMouseMotionHook(GLint x, GLint y);
extern void _glutKeyDownHook(unsigned char key, GLint x, GLint y);
extern void _glutKeyUpHook(unsigned char key, GLint x, GLint y);
extern void _glutSpecialKeyDownHook(GLint key, GLint x, GLint y);
extern void _glutSpecialKeyUpHook(GLint key, GLint x, GLint y);


/*
	*****************************************************************************************************************
	Small and simple vector class. Only thing remotely interesting is the ability to rotate around an arbitrary axis!
	No multiply or divide operations overloaded, as they could have different meanings to what you expect.
	*****************************************************************************************************************
*/
class v3 {
	public:
		float x, y, z;
		
		v3(float _x, float _y, float _z);
		
		v3 operator=(const v3 &v);
		v3 operator+=(v3 &v);
		v3 operator-=(v3 &v);
	
		v3 cross(v3 &v);
		void rotateAround(float theta, v3 axis);
		void normalize();
		void copyTo(float * array);
};


#define JCAMERA_TYPE_FREEMOTION 	0
#define JCAMERA_TYPE_WATCHORIGIN 	1

//	Small and simple camera class. Has a local coordinate frame which is adjusted on orientation rotation etc.
//	Easy to use with gluLookAt() - call see()!
class camera {
	public:
		v3 location, up, look, right;
		int type;
		
		camera();
		virtual ~camera();
		
		/* These operation rotate around local coordinate frame for the camera, ie. around up, look or right. */
		void xrot(float theta);
		void yrot(float theta);
		void zrot(float theta);
		void see();
};



//	Small and simple light class
class light {
	public:
		v3 pos, up, right, look;
		GLfloat light_ambient[4];
		GLfloat light_diffuse[4];
		GLfloat light_specular[4];
		GLint id, show_pos;

		void xrot(float theta);
		void yrot(float theta);
		void zrot(float theta);
		
		light(GLint me);
		void enable();
		void disable();
		void shine();
};


/*
	*****************************************************************************************************************
	GLWindow class. Rather bulky, but it sets everything up for you; it sets the current GLWindow pointer to itself,
	so the glut hooks know where to forward the messages to.
	Simply subclass, and overwrite the non "_glut" functions to get what you want.
	*****************************************************************************************************************
*/
#define GLWINDOW_MAXTIMERS 10
#define GLWINDOW_MAXMOUSEBUTTONS 10
#define GLWINDOW_MAXSPECIALKEYS 256


class GLWindow {
	protected:
		// windowdims = { x, y, width, height }
		unsigned int windowdims[4]; /* x, y, width, height */
		// timers[x][0] = delay, timers[x][1] = reactivate: 1 = yes, 0 = no ( i.e. one-shot timer )
		unsigned int timers[GLWINDOW_MAXTIMERS][2]; 
		unsigned int mousepos[2], mousebuttons[GLWINDOW_MAXMOUSEBUTTONS];	// mbuttons[x] = state.
		unsigned char keys[256], special_keys[GLWINDOW_MAXSPECIALKEYS];		// keys[x] = state.
		
		camera *base_cam;
		light *base_light;
		float factor;
		
		// all data is protected access, so even a single sentinel mutex could be used on every method call.
		bool should_quit_thread;
		pthread_t subthread;
		pthread_attr_t basethread_attr;
		pthread_mutex_t basethread_mutex;
		
	public:
		GLWindow(char *title, GLint x, GLint y, GLint width, GLint height, GLint flags, int argc, char **argv);
		virtual ~GLWindow();
		
		/* methods called by GLUT callback hooks and internal stuff; don't overwrite these unless you know what you're doing! */
		virtual void _glutTimer(GLint value);
		virtual void _glutReshape(GLint width, GLint height);
		virtual void _glutDisplay();
		virtual void _glutMouseButton(GLint button, GLint state, GLint x, GLint y);
		virtual void _glutMouseMotion(GLint x, GLint y);
		
		virtual void _glutKeyDown(unsigned char key, GLint x, GLint y);
		virtual void _glutKeyUp(unsigned char, GLint x, GLint y);

		virtual void _glutSpecialKeyDown(GLint key, GLint x, GLint y);
		virtual void _glutSpecialKeyUp(GLint key, GLint x, GLint y);
		
		/* Non-glut hook methods, which nontheless should not be altered without a good reason! */
		virtual GLint SetTimer(GLint delay, bool repeat);
		virtual int Start(void * (*runme)(void *), int quit_check_every);
		virtual bool ShouldQuitThread();
		virtual void Stop();

		/*	Overwrite these to provide the desired behaviour in subclasses. The functions above pass control straight over to these anyway!
			Note that KeyDown() and KeyUp() could be called by either _glutKeyX() or _glutSpecialKeyX(), so watch out! Separate key state
			arrays are maintained for the key presses of normal or special keys. Override this behaviour if desired. */
		virtual void Timer(GLint value);
		virtual void Reshape(GLint width, GLint height);
		virtual void Display();
		virtual void MouseButton(GLint button, GLint state, GLint x, GLint y);
		virtual void MouseMotion(GLint x, GLint y);
		virtual void KeyDown(GLint key, GLint x, GLint y);
		virtual void KeyUp(GLint key, GLint x, GLint y);
		
		void BaseThreadLock();
		void BaseThreadUnlock();
};


#endif

