#include <QtGui>
#include <QtOpenGL>
#include <QtDebug>

#include "glwidget.h"
#include <QMatrix4x4>

#include "artificialconnections.h"

QString GLWidget::arg(QString argname) {
    int nodespos = qApp->arguments().indexOf(QRegExp("-"+argname+"*"));
    if (nodespos!=-1) {
        return qApp->arguments().at(nodespos+1);
    } else {
        return "";
    }
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    resize(800,800);

    QString filename = qApp->arguments().value(1,"");
    qDebug() << "arglength: " << qApp->arguments().length();
    if (filename == "artificial" || qApp->arguments().length() == 1) {
        qDebug() << "creating new artificial data";
        cons = new ArtificialConnections();

    } else if (arg("nodes")!=""){
        qDebug() << arg("nodes");
        cons = new Connections(arg("nodes"), arg("cons"));
    } else {
        qDebug() << filename;
        cons = new Connections(filename);
    }

    view = new GLfloat[16];
    stuffAlpha = 0.99;
    QVector3D size = cons->max-cons->min;
    float largest = qMax(size.x(),size.y());
    scale = (1/largest)*0.95;
    bg = 1;
    p1 = true;
    p2 = true;

    if (qApp->arguments().indexOf(QRegExp("-writefib"))!=-1) cons->writeBinaryVTK(filename+".fib");
    if (qApp->arguments().indexOf(QRegExp("-writeobj"))!=-1) cons->writeOBJ(filename+".obj");

    if (qApp->arguments().indexOf(QRegExp("-screenshot"))!=-1) screenshot(filename+".png");
    if (qApp->arguments().indexOf(QRegExp("-csv"))!=-1) cons->writeCSVs();
    setFocus();
}

GLWidget::~GLWidget()
{
}

void GLWidget::initializeGL()
{
    selected = new QVector3D(0, 0, 0);
    glEnable(GL_RESCALE_NORMAL);

    glCullFace(GL_FRONT);
    glDisable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);

    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glEnable(GL_LIGHTING);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
    static GLfloat global_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
    static GLfloat specular[] = {0.5f, 0.5f, 0.5f , 1.0f};
    glMateriali(GL_FRONT, GL_SHININESS, 1);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

    static GLfloat lightPosition[4] = { 10000, 10000, 50000, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    glGetFloatv(GL_MODELVIEW_MATRIX, view);

    emit cameraUpdated(view);
    emit cameraUpdatedZoom(scale);

}

void GLWidget::paintGL()
{
    glClearColor(bg,bg,bg,1);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    glScaled(scale,scale,scale);

    glMultMatrixf(view);

    glTranslatef(cons->piv.x(),cons->piv.y(),cons->piv.z());

    glEnable(GL_DEPTH_TEST);

    for (int i = 0; i < cons->edges.size(); ++i) {
        cons->edges[i]->paintGL(intermediateNodes, startAndEndNodes, stuffAlpha, dualGradient, useSpline, cons->selected);
    }

}

void GLWidget::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
    ar = (float)width/(float)height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (ar>1){
        glOrtho(ar*-0.5, ar*+0.5, -0.5, +0.5, -1000, 1000);
    } else {
        glOrtho(-0.5, +0.5, -0.5/ar, +0.5/ar, -1000, 1000);
    }
    glMatrixMode(GL_MODELVIEW);

}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();

    select(event);
}

static bool myUnProject(double winX, double winY, double winZ,
                        const double model[16], const double proj[16], const int viewport[4],
                        double *objX, double *objY, double *objZ)
{
    QMatrix4x4 modelMat, projMat;

    // GL uses column-major, QMatrix4x4 expects row-major, so transpose
    for (int i = 0; i < 16; ++i) {
        modelMat.data()[i] = (float)model[i];
        projMat.data()[i]  = (float)proj[i];
    }

    QMatrix4x4 mvp = projMat * modelMat;
    bool invertible;
    QMatrix4x4 invMvp = mvp.inverted(&invertible);
    if (!invertible)
        return false;

    // normalize window coords to [-1, 1]
    float x = (winX - viewport[0]) / viewport[2] * 2.0f - 1.0f;
    float y = (winY - viewport[1]) / viewport[3] * 2.0f - 1.0f;
    float z = 2.0f * winZ - 1.0f;

    QVector4D in(x, y, z, 1.0f);
    QVector4D out = invMvp * in;

    if (out.w() == 0.0f) return false;

    *objX = out.x() / out.w();
    *objY = out.y() / out.w();
    *objZ = out.z() / out.w();
    return true;
}

bool GLWidget::select(QMouseEvent* event) {
    if (event->button() == Qt::MiddleButton) {
        // Matrices & viewport
        GLdouble modelview[16], projection[16];
        GLint viewport[4];
        GLfloat z;
        GLdouble objx, objy, objz;

        // Use your stored view matrix (float -> double cast)
        for (int i = 0; i < 16; i++)
            modelview[i] = static_cast<GLdouble>(view[i]);

        // Get projection and viewport from OpenGL
        glGetDoublev(GL_PROJECTION_MATRIX, projection);
        glGetIntegerv(GL_VIEWPORT, viewport);

        // Flip y since Qt origin is top-left
        int realY = viewport[3] - event->y();

        // Read depth buffer at mouse position
        glReadPixels(event->x(), realY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);

        // Unproject to object space
        myUnProject(event->x(), realY, z,
                    modelview, projection, viewport,
                    &objx, &objy, &objz);

        qDebug() << "Selected point:" << objx << "," << objy << "," << objz;

        selected->setX(objx);
        selected->setY(objy);
        selected->setZ(objz);
        cons->selectForPoint(selected);

        return true;
    }
    return false;
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    glLoadIdentity();

    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        QMatrix4x4 mat(view[0],view[4],view[8],view[12],view[1],view[5],view[9],view[13],view[2],view[6],view[10],view[14],view[3],view[7],view[11],view[15]);
        QVector3D orig(0, 0, 0);
        QVector3D m = mat.map(orig);
        glTranslatef(m.x(), m.y(), m.z());
        glRotatef(qSqrt(dx*dx+dy*dy)/2.0, dy, dx, 0);
        glTranslatef(-m.x(), -m.y(), -m.z());
    } else if (event->buttons() & Qt::RightButton) {
        glTranslatef(dx/(float)width()*ar/scale, -dy/(float)height()/scale, 0);
    }
    lastPos = event->pos();
    glPushMatrix();
    glMultMatrixf(view);
    glGetFloatv(GL_MODELVIEW_MATRIX, view);
    glPopMatrix();

    emit cameraUpdated(view);
    updateGL();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event){
    cons->sortCons();
    updateGL();
}

void GLWidget::wheelEvent (QWheelEvent *event)
{
    int d = event->angleDelta().y();
    scale *= 1.0-d/1200.0;

    emit cameraUpdatedZoom(scale);

    updateGL();
}

void GLWidget::screenshot(QString name){
    qDebug() << "screenshot";
    QPixmap map = this->renderPixmap(this->width() * 5,this->height() * 5,false);
    QPainter p(&map);
    map.save(name);
}

void GLWidget::save(){
    qDebug() << "start save";
    QString name = "screen.png";
    screenshot(name);
}

void GLWidget::keyPressEvent(QKeyEvent *event) {
    qDebug() << "key:" << event->key();
    if (event->key() == Qt::Key_C) {
        if (bg == 1) {
            bg = 0;
        } else {
            bg =1;
        }
    }
    if (event->key() == Qt::Key_1) p1 = !p1;
    if (event->key() == Qt::Key_2) p2 = !p2;

    if (event->key() == Qt::Key_S) startAndEndNodes = !startAndEndNodes;
    if (event->key() == Qt::Key_I) intermediateNodes = !intermediateNodes;
    if (event->key() == Qt::Key_G) dualGradient = !dualGradient;
    if (event->key() == Qt::Key_U) useSpline = !useSpline;

    updateGL();
}

void GLWidget::stuffSliderChanged(int i){
    stuffAlpha = i/(float)100;
    qDebug() << "stuffAlpha: " << stuffAlpha << Qt::endl;
    updateGL();
}

void GLWidget::updateViewMatrix(GLfloat* matrix) {
    memcpy(view, matrix, sizeof(GLfloat) * 16); // Copy new matrix
    updateGL(); // Request a redraw
}

void GLWidget::setZoomLevel(double value) {
    scale = value;
    updateGL();
}
