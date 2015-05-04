/****************************************************************************
** Meta object code from reading C++ file 'MethodBrowser.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../MethodBrowser.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MethodBrowser.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MethodBrowserCustomButton[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      33,   27,   26,   26, 0x05,

 // slots: signature, parameters, type, tag, flags
      52,   26,   26,   26, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MethodBrowserCustomButton[] = {
    "MethodBrowserCustomButton\0\0index\0"
    "methodClicked(int)\0receiveClick()\0"
};

void MethodBrowserCustomButton::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MethodBrowserCustomButton *_t = static_cast<MethodBrowserCustomButton *>(_o);
        switch (_id) {
        case 0: _t->methodClicked((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->receiveClick(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MethodBrowserCustomButton::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MethodBrowserCustomButton::staticMetaObject = {
    { &QPushButton::staticMetaObject, qt_meta_stringdata_MethodBrowserCustomButton,
      qt_meta_data_MethodBrowserCustomButton, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MethodBrowserCustomButton::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MethodBrowserCustomButton::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MethodBrowserCustomButton::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MethodBrowserCustomButton))
        return static_cast<void*>(const_cast< MethodBrowserCustomButton*>(this));
    return QPushButton::qt_metacast(_clname);
}

int MethodBrowserCustomButton::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void MethodBrowserCustomButton::methodClicked(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_MethodBrowserWidget[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      27,   21,   20,   20, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MethodBrowserWidget[] = {
    "MethodBrowserWidget\0\0index\0"
    "methodClicked(int)\0"
};

void MethodBrowserWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MethodBrowserWidget *_t = static_cast<MethodBrowserWidget *>(_o);
        switch (_id) {
        case 0: _t->methodClicked((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MethodBrowserWidget::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MethodBrowserWidget::staticMetaObject = {
    { &QScrollArea::staticMetaObject, qt_meta_stringdata_MethodBrowserWidget,
      qt_meta_data_MethodBrowserWidget, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MethodBrowserWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MethodBrowserWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MethodBrowserWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MethodBrowserWidget))
        return static_cast<void*>(const_cast< MethodBrowserWidget*>(this));
    return QScrollArea::qt_metacast(_clname);
}

int MethodBrowserWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QScrollArea::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
