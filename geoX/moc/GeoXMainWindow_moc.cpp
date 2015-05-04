/****************************************************************************
** Meta object code from reading C++ file 'GeoXMainWindow.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../GeoXMainWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GeoXMainWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GeoXMainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      16,   15,   15,   15, 0x08,
      42,   15,   15,   15, 0x08,
      70,   15,   15,   15, 0x08,
      96,   15,   15,   15, 0x08,
     127,  109,   15,   15, 0x08,
     157,   15,   15,   15, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_GeoXMainWindow[] = {
    "GeoXMainWindow\0\0on_actionOpen_triggered()\0"
    "on_actionSaveAs_triggered()\0"
    "on_actionExit_triggered()\0initWindow()\0"
    "newExperimentName\0sltExperimentChanged(QString)\0"
    "sltViewerContentChanged()\0"
};

void GeoXMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        GeoXMainWindow *_t = static_cast<GeoXMainWindow *>(_o);
        switch (_id) {
        case 0: _t->on_actionOpen_triggered(); break;
        case 1: _t->on_actionSaveAs_triggered(); break;
        case 2: _t->on_actionExit_triggered(); break;
        case 3: _t->initWindow(); break;
        case 4: _t->sltExperimentChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->sltViewerContentChanged(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData GeoXMainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject GeoXMainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_GeoXMainWindow,
      qt_meta_data_GeoXMainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &GeoXMainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *GeoXMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *GeoXMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GeoXMainWindow))
        return static_cast<void*>(const_cast< GeoXMainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int GeoXMainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
