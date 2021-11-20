#include "gl.h"

GLWindow * _currentGLWindow = NULL;


// VECTOR CLASS
v3::v3(float _x, float _y, float _z) { x = _x; y = _y; z = _z; }
v3 v3::operator=(const v3 &v) { if (this==&v) return *this; x = v.x; y = v.y; z = v.z; return *this; }
v3 v3::operator+=(v3 &v) { x += v.x; y += v.y; z += v.z; return *this; }
v3 v3::operator-=(v3 &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
v3 v3::cross(v3 &v) { return v3( y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x ); }


void v3::rotateAround(float theta, v3 axis)
{
	float nx = 0.0, ny = 0.0, nz = 0.0;
	float cos_theta = cos(theta), sin_theta = sin(theta);
	
	/* assumes axis is normalized! */
	nx += (cos_theta + (1 - cos_theta) * axis.x * axis.x) * x;
	nx += ((1 - cos_theta) * axis.x * axis.y - axis.z * sin_theta) * y;
	nx += ((1 - cos_theta) * axis.x * axis.z + axis.y * sin_theta) * z;

	ny += ((1 - cos_theta) * axis.x * axis.y + axis.z * sin_theta) * x;
	ny += (cos_theta + (1 - cos_theta) * axis.y * axis.y) * y;
	ny += ((1 - cos_theta) * axis.y * axis.z - axis.x * sin_theta) * z;

	nz += ((1 - cos_theta) * axis.x * axis.z - axis.y * sin_theta) * x;
	nz += ((1 - cos_theta) * axis.y * axis.z + axis.x * sin_theta) * y;
	nz += (cos_theta + (1 - cos_theta) * axis.z * axis.z) * z;
	
	x = nx;
	y = ny;
	z = nz;
}


void v3::normalize()
{
	float m = sqrt( x*x + y*y + z*z );
	x /= m;
	y /= m;
	z /= m;
}


void v3::copyTo( float * array )
{
	if( array != NULL )
	{
		array[0] = x;
		array[1] = y;
		array[2] = z;
	}
}





/*
	*****************************************************************************************************************
	Small and simple camera class. Has a local coordinate frame which is adjusted on orientation rotation etc.
	Easy to use with gluLookAt() - call see()!
	*****************************************************************************************************************
*/
camera::camera() : location(0.0, 0.0, 0.0), up(0.0, 1.0, 0.0), look(0.0, 0.0, 1.0), right(1.0, 0.0, 0.0) { type = JCAMERA_TYPE_FREEMOTION; }


camera::~camera() {}

void camera::xrot(float theta) {
	if (type == JCAMERA_TYPE_FREEMOTION) {
		/*	If we assume the axis of rotation is normalized, then we only have to rotate one of the vectors
			as we can recalculate the other from the cross product of the rotated vector and the axis of rotation.
			this removes costly tri functions, and possibly a sqrt() from a normalize call */
		up.rotateAround( theta, right );
		up.normalize();
		look = right.cross( up );
	}
	else if (type == JCAMERA_TYPE_WATCHORIGIN) {
		location.rotateAround( theta, right );
		look.x = -location.x;
		look.y = -location.y;
		look.z = -location.z;
		look.normalize();
		up = look.cross(right);
	}
}


void camera::yrot(float theta) {
	if (type == JCAMERA_TYPE_FREEMOTION) {
		right.rotateAround(theta, up);
		right.normalize();
		look = right.cross(up);
	}
	else if (type == JCAMERA_TYPE_WATCHORIGIN) {
		location.rotateAround(theta, up);
		look.x = -location.x;
		look.y = -location.y;
		look.z = -location.z;
		look.normalize();
		right = up.cross(look);
	}
}


void camera::zrot(float theta) {
	// do same thing regardless of type!
	up.rotateAround(theta, look);
	up.normalize();
	right = up.cross(look);
}


void camera::see()
{
	gluLookAt(location.x, location.y, location.z, 
		       location.x + look.x, location.y + look.y, location.z + look.z, 
			   up.x, up.y, up.z);
}




/*
	*****************************************************************************************************************
	Small and simple light class.
	*****************************************************************************************************************
*/

light::light(GLint me) : pos(0.0, 0.0, 0.0), up(0.0, 1.0, 0.0), right(1.0, 0.0, 0.0), look(0.0, 0.0, 1.0) {

	id = me;
	show_pos = true;
	
	light_ambient[0] = 0.0;
	light_ambient[1] = 0.0;
	light_ambient[2] = 0.0;
	light_ambient[3] = 1.0;

	light_diffuse[0] = 1.0;
	light_diffuse[1] = 1.0;
	light_diffuse[2] = 1.0;
	light_diffuse[3] = 1.0;

	light_specular[0] = 1.0;
	light_specular[1] = 1.0;
	light_specular[2] = 1.0;
	light_specular[3] = 1.0;
}


void light::xrot(float theta) {
	pos.rotateAround(theta, right);
	look.x = -pos.x;
	look.y = -pos.y;
	look.z = -pos.z;
	look.normalize();
	up = look.cross(right);
}


void light::yrot(float theta) {
	pos.rotateAround(theta, up);
	look.x = -pos.x;
	look.y = -pos.y;
	look.z = -pos.z;
	look.normalize();
	right = up.cross(look);
}


void light::zrot(float theta) {
	up.rotateAround(theta, look);
	up.normalize();
	right = up.cross(look);
}


void light::enable() { glEnable(id); }


void light::disable()	{ glDisable(id); }


void light::shine() {
	GLfloat light_pos[4];

	glLightfv( id, GL_AMBIENT,  light_ambient );
	glLightfv( id, GL_DIFFUSE,  light_diffuse );
	glLightfv( id, GL_SPECULAR, light_specular );
	
	pos.copyTo( light_pos );
	light_pos[3] = 1.0;
	
	glPushMatrix();
	glLightfv( id, GL_POSITION, light_pos );
	glPopMatrix();
}



/*
	*****************************************************************************************************************
	GLWindow class. Rather bulky, but it sets everything up for you; it sets the current GLWindow pointer to itself,
	so the glut hooks know where to forward the messages to.
	Simply subclass, and overwrite the non "_glut" functions to get what you want.
	*****************************************************************************************************************
*/
GLWindow::GLWindow(char *title, GLint x, GLint y, GLint width, GLint height, GLint flags, int argc, char ** argv) {

	int i;
	
	windowdims[0] = x;
	windowdims[1] = y;
	windowdims[2] = width;
	windowdims[3] = height;
		
	SetCurrentGLWindow(this);

	glutInit(&argc, argv);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(flags);
	glutCreateWindow(title);
		
	glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);

	/*	Set up hook functions; these pass the message on to the current GLWindow instance _glutX methods, which in turn pass
		the messages onto the appropriate virtual methods, overwritten by subclasses. Could do this outside any instantiation,
		but do it automatically here.
		NOTE: would these cause problems with threaded code and muliple GLWindow instances? Re-entrant? */
	glutDisplayFunc( _glutDisplayHook );
	glutReshapeFunc( _glutReshapeHook );
	glutMouseFunc( _glutMouseButtonHook );
	glutMotionFunc( _glutMouseMotionHook );
	glutKeyboardFunc( _glutKeyDownHook );
	glutKeyboardUpFunc( _glutKeyUpHook );
	glutSpecialFunc( _glutSpecialKeyDownHook );
	glutSpecialUpFunc( _glutSpecialKeyUpHook );

	base_cam = new camera();
	base_light = new light( GL_LIGHT0 );
	base_light->enable();
	
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	
	// wipe timer info - could cause weirdness otherwise
	for (i = 0; i < GLWINDOW_MAXTIMERS; i++) {
		timers[i][0] = 0;
		timers[i][1] = 0;
	}
}


GLWindow::~GLWindow() { delete base_cam; delete base_light; }


void GLWindow::_glutTimer( GLint value )
{
	int rv;
	void * vptr;
	
	/*
		Is this the quit timer? Should we stop here? Attempt to join the child thread.
	*/
	if( value == 0 && ShouldQuitThread() == true)
	{
//		printf( "%s(): stop timer recieved - processing...\n", __func__ );
		/* allow cleanup in child thread before calling derived Timer(), which might rely on the child thread in some way. */
		rv = pthread_join( subthread, &vptr );
		if (rv)
		{
			if( rv == EINVAL ) printf( "%s(): ERROR; return code from pthread_join() is %d ( EINVAL )\n", __func__, rv );
			else if( rv == ESRCH ) printf( "%s(): ERROR; return code from pthread_join() is %d ( ESRCH )\n", __func__, rv );
			else if( rv == EDEADLK ) printf( "%s(): ERROR; return code from pthread_join() is %d ( EDEADLK )\n", __func__, rv );
		}
//		else printf("%s(): completed join with thread status = %d\n", __func__, rv);
		Timer(value);
		delete this;
		exit(0);
	}
	else
	{
		Timer(value);
		/* automatically resart timer if it's set to be repeating. */
		if( timers[value][1] == 1 ) glutTimerFunc( timers[value][0], _glutTimerHook, value );
		/* otherwise mark as free, so we can recycle later in SetTimer() */
		else
		{
			BaseThreadLock();
			timers[value][0] = 0;
			timers[value][1] = 0;
			BaseThreadUnlock();
		}
	}
}


int GLWindow::Start(void * (*runme)(void *), int quit_check_every)
{
	/*
		This runs the runme() funtion in a child thread. It would be much better to have the actual visualisation window as a child
		thread, but Apple and others are muppets and only allow GUI messages to be delivered to the main thread, so I have to jump
		through hoops instead of having a SENSIBLE heirarchy.
	*/
	int rv;
	if( (rv = pthread_mutex_init( &basethread_mutex, NULL )) ) { printf( "%s(): unable to make mutex ( error code %d ).\n", __func__, rv ); delete this; exit(-1); }
	if( (rv = pthread_attr_init( &basethread_attr )) )  { printf( "%s(): unable to make pthread attr ( error code %d ).\n", __func__, rv ); delete this; exit(-1); }
	if( (rv = pthread_attr_setdetachstate( &basethread_attr, PTHREAD_CREATE_JOINABLE )) ) { printf( "%s(): unable to make pthread joinable ( error code %d ).\n", __func__, rv ); delete this; exit(-1); }
	if( (rv = pthread_create( &subthread, &basethread_attr, runme, this )) ) { printf( "%s(): unable to make pthread ( error code %d ).\n", __func__, rv ); delete this; exit(-1); }
	/* Quit timer */
	should_quit_thread = false;
	timers[0][0] = quit_check_every;
	timers[0][1] = 1;
	glutTimerFunc( quit_check_every, _glutTimerHook, 0 );
	glutMainLoop();
	return rv;
}


bool GLWindow::ShouldQuitThread()
{
	return should_quit_thread;
}


void GLWindow::Stop()
{
//	printf( "%s(): recieved stop message.\n", __func__ );
	BaseThreadLock();
	should_quit_thread = true;
	BaseThreadUnlock();
}


void GLWindow::_glutReshape(GLint width, GLint height) {
	windowdims[2] = width;
	windowdims[3] = height;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, width, height);
	if (ORTHOGRAPHIC) glOrtho(-width*factor, width*factor, -height*factor, height*factor, 0.01, 1000000.0);
	else gluPerspective(45, (float) width / height, 1, 100000);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	Reshape(width, height);
}


void GLWindow::_glutDisplay() {
	Display();
}


void GLWindow::_glutMouseButton(GLint button, GLint state, GLint x, GLint y) {
	if (button >= GLWINDOW_MAXMOUSEBUTTONS) {
		printf( "Mouse button %d pressed; GLWINDOW_MAXMOUSEBUTTONS is %d, so increase! This press is ignored.\n", button, GLWINDOW_MAXMOUSEBUTTONS );
		return;
	}
	mousebuttons[button] = state;
	mousepos[0] = x;
	mousepos[1] = y;
	MouseButton(button, state, x, y);
}


void GLWindow::_glutMouseMotion(GLint x, GLint y) {
	mousepos[0] = x;
	mousepos[1] = y;
	MouseMotion(x, y);
}


void GLWindow::_glutKeyDown(unsigned char key, GLint x, GLint y) {
	keys[key] = 1;
	KeyDown((GLint)key, x, y);
}


void GLWindow::_glutKeyUp(unsigned char key, GLint x, GLint y) {
	keys[key] = 0;
	KeyUp((GLint)key, x, y);
}


void GLWindow::_glutSpecialKeyDown(GLint key, GLint x, GLint y) {
	if (key > GLWINDOW_MAXSPECIALKEYS-1) {
		printf("Key %d (%c) pressed; GLWINDOW_MAXSPECIALKEYS is %d, so increase! This press is ignored.\n", key, key, GLWINDOW_MAXSPECIALKEYS);
		return;
	}
	special_keys[key] = 1;
	KeyDown(key, x, y);
}


void GLWindow::_glutSpecialKeyUp(GLint key, GLint x, GLint y) {
	if(key > GLWINDOW_MAXSPECIALKEYS-1) {
		printf("Key %d (%c) released; GLWINDOW_MAXSPECIALKEYS is %d, so increase! This release is ignored.\n", key, key, GLWINDOW_MAXSPECIALKEYS);
		return;
	}
	special_keys[key] = 0;
	KeyUp(key, x, y);
}


GLint GLWindow::SetTimer(GLint delay, bool repeat) {
	int i;
	/* zero-delay timers not allowed; why not just call the function directly? */
	if (delay == 0) return -1;
	/*
		Find first empty slot in timer array.
		Note that the timer id doubles as the index into the timer array; hence, always unique id reurned where success.
		There are probably more efficient ways to do this rather than iteraing over the array, but this should be fast enough.
		
		index 0 is special timer for quit notification, so start at 1.
	*/
	for (i=1; i<GLWINDOW_MAXTIMERS; i++) {
		if ((timers[i][0] == 0) && (timers[i][1] == 0)) {
			BaseThreadLock();
			timers[i][0] = delay;
			timers[i][1] = (repeat == true ? 1 : 0);
			BaseThreadUnlock();
			glutTimerFunc(delay, _glutTimerHook, i);
			return i;
		}
	}
	printf( "Unable to create timer! Is GLWINDOW_MAXTIMER too small, or hasn't be initialised as zeros?\n" );
	/* unable to fulfil request for timer! */
	return -1;
}

void GLWindow::Timer(GLint value) {}
void GLWindow::Reshape(GLint width, GLint height) { Display(); }
void GLWindow::Display() {}
void GLWindow::MouseButton(GLint button, GLint state, GLint x, GLint y) {}
void GLWindow::MouseMotion(GLint x, GLint y) {}
void GLWindow::KeyDown(GLint key, GLint x, GLint y) {}
void GLWindow::KeyUp(GLint key, GLint x, GLint y) {}

void GLWindow::BaseThreadLock() {  pthread_mutex_lock(&basethread_mutex);  }
void GLWindow::BaseThreadUnlock() {  pthread_mutex_unlock(&basethread_mutex);  }


/*
	*****************************************************************************************************************
	These are the c-callable functions that are passed to the glut callback routines. Don't bother with these,
	they're boring - all they do is check for a valid current GLWindow to forward messages to!
	There's also an error function in here, and a function to set the current GLWindow pointer.
	*****************************************************************************************************************
*/

void ErrorMessage(const char *calledfrom, const int linenumber, const char *msg) {
	printf("%s ( line %d ): %s\n", calledfrom, linenumber, msg);
	exit(-1);
}

void SetCurrentGLWindow(GLWindow *glw) {
	if( glw == NULL ) ErrorMessage(__FILE__, __LINE__, "Null pointer to GLWindow passed.");
	_currentGLWindow = glw;
}

void _glutTimerHook(GLint value) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutTimer(value);
}

void _glutReshapeHook(GLint width, GLint height) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutReshape(width, height );
}

void _glutDisplayHook() {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutDisplay();
}

void _glutMouseButtonHook(GLint button, GLint state, GLint x, GLint y) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutMouseButton(button, state, x, y);
}

void _glutMouseMotionHook(GLint x, GLint y) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutMouseMotion(x, y);
}

void _glutKeyDownHook(unsigned char key, GLint x, GLint y) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutKeyDown((GLint)key, x, y);
}

void _glutKeyUpHook(unsigned char key, GLint x, GLint y) {
	if( _currentGLWindow == NULL ) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutKeyUp((GLint)key, x, y);
}

void _glutSpecialKeyDownHook(GLint key, GLint x, GLint y) {
	if(_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutSpecialKeyDown(key, x, y);
}

void _glutSpecialKeyUpHook(GLint key, GLint x, GLint y) {
	if (_currentGLWindow == NULL) ErrorMessage(__FILE__, __LINE__, "No acive GLWindow.");
	_currentGLWindow->_glutSpecialKeyUp(key, x, y);
}

