/****************************************************************************
** Meta object code from reading C++ file 'GeoXGLWidget3D.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../GeoXGLWidget3D.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GeoXGLWidget3D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GeoXGLWidget3D[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      16,   15,   15,   15, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_GeoXGLWidget3D[] = {
    "GeoXGLWidget3D\0\0widgetContentChanged()\0"
};

void GeoXGLWidget3D::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        GeoXGLWidget3D *_t = static_cast<GeoXGLWidget3D *>(_o);
        switch (_id) {
        case 0: _t->widgetContentChanged(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData GeoXGLWidget3D::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject GeoXGLWidget3D::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GeoXGLWidget3D,
      qt_meta_data_GeoXGLWidget3D, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &GeoXGLWidget3D::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *GeoXGLWidget3D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *GeoXGLWidget3D::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GeoXGLWidget3D))
        return static_cast<void*>(const_cast< GeoXGLWidget3D*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int GeoXGLWidget3D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void GeoXGLWidget3D::widgetContentChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
