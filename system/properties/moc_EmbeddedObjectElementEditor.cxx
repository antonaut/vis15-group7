/****************************************************************************
** Meta object code from reading C++ file 'EmbeddedObjectElementEditor.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "EmbeddedObjectElementEditor.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EmbeddedObjectElementEditor.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_EmbeddedObjectElementEditor[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      29,   28,   28,   28, 0x0a,
      45,   41,   28,   28, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_EmbeddedObjectElementEditor[] = {
    "EmbeddedObjectElementEditor\0\0sltButton()\0"
    "obj\0sltDestroyEditorWidget(QObject*)\0"
};

void EmbeddedObjectElementEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        EmbeddedObjectElementEditor *_t = static_cast<EmbeddedObjectElementEditor *>(_o);
        switch (_id) {
        case 0: _t->sltButton(); break;
        case 1: _t->sltDestroyEditorWidget((*reinterpret_cast< QObject*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData EmbeddedObjectElementEditor::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject EmbeddedObjectElementEditor::staticMetaObject = {
    { &ElementEditor::staticMetaObject, qt_meta_stringdata_EmbeddedObjectElementEditor,
      qt_meta_data_EmbeddedObjectElementEditor, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &EmbeddedObjectElementEditor::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *EmbeddedObjectElementEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *EmbeddedObjectElementEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_EmbeddedObjectElementEditor))
        return static_cast<void*>(const_cast< EmbeddedObjectElementEditor*>(this));
    return ElementEditor::qt_metacast(_clname);
}

int EmbeddedObjectElementEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = ElementEditor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}
static const uint qt_meta_data_EmbeddedObjectInlineElementEditor[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

static const char qt_meta_stringdata_EmbeddedObjectInlineElementEditor[] = {
    "EmbeddedObjectInlineElementEditor\0"
};

void EmbeddedObjectInlineElementEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObjectExtraData EmbeddedObjectInlineElementEditor::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject EmbeddedObjectInlineElementEditor::staticMetaObject = {
    { &ElementEditor::staticMetaObject, qt_meta_stringdata_EmbeddedObjectInlineElementEditor,
      qt_meta_data_EmbeddedObjectInlineElementEditor, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &EmbeddedObjectInlineElementEditor::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *EmbeddedObjectInlineElementEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *EmbeddedObjectInlineElementEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_EmbeddedObjectInlineElementEditor))
        return static_cast<void*>(const_cast< EmbeddedObjectInlineElementEditor*>(this));
    return ElementEditor::qt_metacast(_clname);
}

int EmbeddedObjectInlineElementEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = ElementEditor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
QT_END_MOC_NAMESPACE
