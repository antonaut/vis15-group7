/****************************************************************************
** Meta object code from reading C++ file 'ClassTypeClassPropertyEditorWidget.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ClassTypeClassPropertyEditorWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ClassTypeClassPropertyEditorWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ClassTypeClassPropertyEditorWidget[] = {

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
      36,   35,   35,   35, 0x05,

 // slots: signature, parameters, type, tag, flags
      56,   50,   35,   35, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ClassTypeClassPropertyEditorWidget[] = {
    "ClassTypeClassPropertyEditorWidget\0\0"
    "sigModified()\0index\0selectionChanged(int)\0"
};

void ClassTypeClassPropertyEditorWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ClassTypeClassPropertyEditorWidget *_t = static_cast<ClassTypeClassPropertyEditorWidget *>(_o);
        switch (_id) {
        case 0: _t->sigModified(); break;
        case 1: _t->selectionChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ClassTypeClassPropertyEditorWidget::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ClassTypeClassPropertyEditorWidget::staticMetaObject = {
    { &QComboBox::staticMetaObject, qt_meta_stringdata_ClassTypeClassPropertyEditorWidget,
      qt_meta_data_ClassTypeClassPropertyEditorWidget, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ClassTypeClassPropertyEditorWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ClassTypeClassPropertyEditorWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ClassTypeClassPropertyEditorWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ClassTypeClassPropertyEditorWidget))
        return static_cast<void*>(const_cast< ClassTypeClassPropertyEditorWidget*>(this));
    return QComboBox::qt_metacast(_clname);
}

int ClassTypeClassPropertyEditorWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QComboBox::qt_metacall(_c, _id, _a);
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
void ClassTypeClassPropertyEditorWidget::sigModified()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
