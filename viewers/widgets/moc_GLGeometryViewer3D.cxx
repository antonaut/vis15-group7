/****************************************************************************
** Meta object code from reading C++ file 'GLGeometryViewer3D.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GLGeometryViewer3D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GLGeometryViewer3D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GLGeometryViewer3D[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      20,   19,   19,   19, 0x05,

 // slots: signature, parameters, type, tag, flags
      43,   19,   19,   19, 0x08,
      66,   19,   19,   19, 0x08,
      87,   19,   19,   19, 0x08,
     109,   19,   19,   19, 0x08,
     131,   19,   19,   19, 0x08,
     156,   19,   19,   19, 0x08,
     181,   19,   19,   19, 0x08,
     205,   19,   19,   19, 0x08,
     233,   19,   19,   19, 0x08,
     264,   19,   19,   19, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_GLGeometryViewer3D[] = {
    "GLGeometryViewer3D\0\0viewerContentChanged()\0"
    "on_btnCamera_pressed()\0on_btnPick_pressed()\0"
    "on_btnClear_pressed()\0on_btnLight_pressed()\0"
    "on_btnFillMode_pressed()\0"
    "on_btnDrawAxis_pressed()\0"
    "on_btnHandles_pressed()\0"
    "on_btnResetCamera_pressed()\0"
    "on_btnVertexLighting_pressed()\0"
    "sltWidgetContentChanged()\0"
};

void GLGeometryViewer3D::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        GLGeometryViewer3D *_t = static_cast<GLGeometryViewer3D *>(_o);
        switch (_id) {
        case 0: _t->viewerContentChanged(); break;
        case 1: _t->on_btnCamera_pressed(); break;
        case 2: _t->on_btnPick_pressed(); break;
        case 3: _t->on_btnClear_pressed(); break;
        case 4: _t->on_btnLight_pressed(); break;
        case 5: _t->on_btnFillMode_pressed(); break;
        case 6: _t->on_btnDrawAxis_pressed(); break;
        case 7: _t->on_btnHandles_pressed(); break;
        case 8: _t->on_btnResetCamera_pressed(); break;
        case 9: _t->on_btnVertexLighting_pressed(); break;
        case 10: _t->sltWidgetContentChanged(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData GLGeometryViewer3D::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject GLGeometryViewer3D::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_GLGeometryViewer3D,
      qt_meta_data_GLGeometryViewer3D, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &GLGeometryViewer3D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *GLGeometryViewer3D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *GLGeometryViewer3D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GLGeometryViewer3D))
        return static_cast<void*>(const_cast< GLGeometryViewer3D*>(this));
    return QWidget::qt_metacast(_clname);
}

int GLGeometryViewer3D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    }
    return _id;
}

// SIGNAL 0
void GLGeometryViewer3D::viewerContentChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
